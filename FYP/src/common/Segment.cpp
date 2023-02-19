#include "Segment.h"
#include <string>
#include "utils.h"


FYP::Segment::Segment(unsigned int original_id_, Point p1, Point p2) {
	this->original_id_ = original_id_;
	this->start = p1;
	this->end = p2;
}

void FYP::Segment::setSegment(Point& p1, Point& p2) {
	memcpy(start.getDataAddr(), p1.getDataAddr(), sizeof(double) * 3);
	memcpy(end.getDataAddr(), p2.getDataAddr(), sizeof(double) * 3);
}

unsigned int FYP::Segment::serialize(char* content) {
	unsigned int pos = 0;
	memcpy(content + pos, &original_id_, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(content + pos, start.getDataAddr(), sizeof(double) * 3);
	pos += sizeof(double) * 3;
	memcpy(content + pos, end.getDataAddr(), sizeof(double) * 3);
	pos += sizeof(double) * 3;
	return pos;
}

unsigned int FYP::Segment::deserialize(char* content) {
	unsigned int pos = 0;
	memcpy(&original_id_, content + pos, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(start.getDataAddr(), content + pos, sizeof(double) * 3);
	pos += sizeof(double) * 3;
	memcpy(end.getDataAddr(), content + pos, sizeof(double) * 3);
	pos += sizeof(double) * 3;
	return pos;
}

unsigned int FYP::Segment::getByteSize() const {
	unsigned int size = 0;
	size += sizeof(original_id_);
	size += sizeof(double) * 6;
	return size;
}

bool FYP::Segment::isOverlap(double query_max[3], double query_min[3]) {
	bool is_overlap1 = isLineRectangleOverlap(start.getX(), start.getY(), end.getX(), end.getY(), query_min[0], query_min[1], query_max[0], query_max[1]);
	bool is_overlap2 = isLineRectangleOverlap(start.getX(), start.getT(), end.getX(), end.getT(), query_min[0], query_min[2], query_max[0], query_max[2]);
	bool is_overlap3 = isLineRectangleOverlap(start.getY(), start.getT(), end.getY(), end.getT(), query_min[1], query_min[2], query_max[1], query_max[2]);
	bool is_overlap = is_overlap1 && is_overlap2 && is_overlap3;
	return is_overlap;
}

bool FYP::Segment::isTemporalOverlap(double temporal_range[2]) {
	if (start.getT() < temporal_range[0] || end.getT() > temporal_range[1]) {
		return false;
	}
	return true;
}
