#pragma once

#include "Point.h"

namespace FYP {
	class Segment {
	private:
		unsigned int original_id_;
		Point start;
		Point end;
	public:
		Segment(){};
		Segment(unsigned int original_id_, Point p1, Point p2);
		void setSegment(Point& p1, Point& p2);
		Point* getStart() { return &start; }
		Point* getEnd() { return &end; }
		unsigned int* getId() { return &original_id_; }
		unsigned int serialize(char* content);
		unsigned int deserialize(char* content);
		unsigned int getByteSize() const;
		unsigned int getOriginalId() const { return original_id_; }
		void setOriginalId(unsigned int original_id_) { this->original_id_ = original_id_; }
		bool isOverlap(double query_max[3], double query_min[3]);
		bool isTemporalOverlap(double temporal_range[2]);
	};
}

