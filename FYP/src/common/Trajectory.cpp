#include "Trajectory.h"
#include "utils.h"
#include <math.h>
#include <vector>
#include <iostream>

FYP::Trajectory::Trajectory(unsigned int id, unsigned int point_num) {
	original_id_ = id;
	id_ = id;
	points_ = std::vector<Point>(point_num);
}

unsigned int FYP::Trajectory::getPageSize() const {
	unsigned int size = sizeof(double) * 3 * points_.size() + sizeof(unsigned int);
	return size;
}

void FYP::Trajectory::addPoint(const double x, const double y, const double t) {
	points_.emplace_back(x, y, t);
}

void FYP::Trajectory::addPoint(const Point& point) {
	points_.emplace_back(point.getX(), point.getY(), point.getT());
}

bool FYP::Trajectory::isOverlap(const unsigned int dim, const double value) const {
	for (int i = 0; i < points_.size() - 1; i++) {
		double max_value = std::max(points_[i].getValueFromDim(dim), points_[i + 1].getValueFromDim(dim));
		double min_value = std::min(points_[i].getValueFromDim(dim), points_[i + 1].getValueFromDim(dim));
		if (min_value <= value && value < max_value) {
			return true;
		}
	}
	return false;
}

bool FYP::Trajectory::isOverlap(double* query_max, double* query_min) {
	for (unsigned int i = 0; i < points_.size() - 1; i++) {
		bool is_overlap = checkOverlap(points_[i], points_[i + 1], query_max, query_min);
		if (is_overlap) {
			return true;
		}
	}
	return false;
}

bool FYP::Trajectory::isTemporalOverlap(double temporal_range[2]) {
	if (points_.front().getT() < temporal_range[0] || points_.back().getT() > temporal_range[1]) {
		return false;
	}
	return true;
}

void FYP::Trajectory::crossSplit(const double cut[2], std::vector<std::list<Trajectory*> >& sub_trajectories) {
	//cout << "Cross splitting trajectory with" << cut[0]<<" "<<cut[1]<<endl;
	//Print();
	int start_pos[2] = { -1, -1 };
	if (points_[0].getX() < cut[0]) {
		start_pos[0] = 0;
	}
	if (points_[0].getX() > cut[0]) {
		start_pos[0] = 1;
	}
	if (points_[0].getY() < cut[1]) {
		start_pos[1] = 0;
	}
	if (points_[0].getY() > cut[1]) {
		start_pos[1] = 1;
	}
	//cout << "start_pos: " << start_pos[0] << " " << start_pos[1] << endl;
	int sub_traj_id = 0;
	int next_pos[2] = { -1, -1 };
	Trajectory* sub_traj = new Trajectory();
	sub_traj->id_ = sub_traj_id;
	sub_traj->original_id_ = original_id_;
	sub_traj->points_.push_back(points_[0]);
	int sub_area_id = 0;
	for (unsigned int i = 1; i < points_.size(); i++) {
		if (points_[i].getX() < cut[0]) {
			next_pos[0] = 0;
		}
		else if (points_[i].getX() > cut[0]) {
			next_pos[0] = 1;
		}
		else {
			next_pos[0] = -1;
		}

		if (points_[i].getY() < cut[1]) {
			next_pos[1] = 0;
		}
		else if (points_[i].getY() > cut[1]) {
			next_pos[1] = 1;
		}
		else {
			next_pos[1] = -1;
		}
		//cout << "next_pos: " << next_pos[0] << " " << next_pos[1] << endl;
		if (start_pos[0] + next_pos[0] != 1 && start_pos[1] + next_pos[1] != 1) {
			sub_traj->points_.push_back(points_[i]);
			if (start_pos[0] == -1) {
				start_pos[0] = next_pos[0];
			}
			if (start_pos[1] == -1) {
				start_pos[1] = next_pos[1];
			}
		}
		else if (start_pos[0] + next_pos[0] == 1 && start_pos[1] + next_pos[1] == 1) {
			double ratio1 = (cut[0] - points_[i - 1].getX()) / (points_[i].getX() - points_[i - 1].getX());
			double ratio2 = (cut[1] - points_[i - 1].getY()) / (points_[i].getY() - points_[i - 1].getY());
			double small_ratio = std::min(ratio1, ratio2);
			double large_ratio = std::max(ratio1, ratio2);
			//cout << "ratio: " << small_ratio << " " << large_ratio << endl;
			//cout << (small_ratio == 0) << " " << (large_ratio == 0) << endl;
			double coord1[3]{};
			double coord2[3]{};
			coord1[0] = points_[i - 1].getX() + small_ratio * (points_[i].getX() - points_[i - 1].getX());
			coord1[1] = points_[i - 1].getY() + small_ratio * (points_[i].getY() - points_[i - 1].getY());
			coord1[2] = points_[i - 1].getT() + small_ratio * (points_[i].getT() - points_[i - 1].getT());
			coord2[0] = points_[i - 1].getX() + large_ratio * (points_[i].getX() - points_[i - 1].getX());
			coord2[1] = points_[i - 1].getY() + large_ratio * (points_[i].getY() - points_[i - 1].getY());
			coord2[2] = points_[i - 1].getT() + large_ratio * (points_[i].getT() - points_[i - 1].getT());
			Point p1(coord1);
			Point p2(coord2);
			sub_area_id = start_pos[0] + start_pos[1] * 2;
			if (small_ratio == 0) {
				sub_trajectories[sub_area_id].push_back(sub_traj);
				sub_traj_id += 1;
				sub_traj = new Trajectory();
				sub_traj->id_ = sub_traj_id;
				sub_traj->original_id_ = original_id_;
				sub_traj->points_.push_back(points_[i - 1]);
				if (large_ratio == 0) {
					sub_traj->points_.push_back(points_[i]);
				}
				else {
					sub_traj->points_.push_back(p2);
					sub_area_id = 0;
					if (0.5 * (p1.getX() + p2.getX()) > cut[0]) {
						sub_area_id = 1;
					}
					if (0.5 * (p1.getY() + p2.getY()) > cut[1]) {
						sub_area_id += 2;
					}
					sub_trajectories[sub_area_id].push_back(sub_traj);
					sub_traj_id += 1;
					sub_traj = new Trajectory();
					sub_traj->id_ = sub_traj_id;
					sub_traj->original_id_ = original_id_;
					sub_traj->points_.push_back(p2);
					sub_traj->points_.push_back(points_[i]);
				}

			}
			else {
				sub_traj->points_.push_back(p1);
				sub_area_id = start_pos[0] + start_pos[1] * 2;
				sub_trajectories[sub_area_id].push_back(sub_traj);
				sub_traj_id += 1;
				sub_traj = new Trajectory();
				sub_traj->original_id_ = original_id_;
				sub_traj->id_ = sub_traj_id;
				sub_traj->points_.push_back(p1);
				sub_traj->points_.push_back(p2);
				sub_area_id = 0;
				if (0.5 * (p1.getX() + p2.getX()) > cut[0]) {
					sub_area_id = 1;
				}
				if (0.5 * (p1.getY() + p2.getY()) > cut[1]) {
					sub_area_id += 2;
				}
				sub_trajectories[sub_area_id].push_back(sub_traj);
				sub_traj_id += 1;
				sub_traj = new Trajectory();
				sub_traj->id_ = sub_traj_id;
				sub_traj->original_id_ = original_id_;
				sub_traj->points_.push_back(p2);
				sub_traj->points_.push_back(points_[i]);
			}
			start_pos[0] = next_pos[0];
			start_pos[1] = next_pos[1];

		}
		else {
			int dim = 0;
			if (start_pos[1] + next_pos[1] == 1) {
				dim = 1;
			}
			double ratio = (cut[dim] - points_[i - 1].getValueFromDim(dim)) / (points_[i].getValueFromDim(dim) - points_[i - 1].getValueFromDim(dim));
			//cout << "ratio: " << ratio << endl;
			double coord[3];
			coord[0] = points_[i - 1].getX() + ratio * (points_[i].getX() - points_[i - 1].getX());
			coord[1] = points_[i - 1].getY() + ratio * (points_[i].getY() - points_[i - 1].getY());
			coord[2] = points_[i - 1].getT() + ratio * (points_[i].getT() - points_[i - 1].getT());
			Point p(coord);
			sub_area_id = 0;
			if (start_pos[0] == 1) {
				sub_area_id = 1;
			}
			if (start_pos[1] == 1) {
				sub_area_id += 2;
			}
			if (ratio == 0) {
				//sub_traj->points_.push_back(points_[i]);
			}
			else {

				sub_traj->points_.push_back(p);
			}
			sub_trajectories[sub_area_id].push_back(sub_traj);
			sub_traj_id += 1;
			sub_traj = new Trajectory();
			sub_traj->id_ = sub_traj_id;
			sub_traj->original_id_ = original_id_;
			if (ratio == 0) {
				sub_traj->points_.push_back(points_[i - 1]);
				sub_traj->points_.push_back(points_[i]);
				if (points_[i - 1].getX() < cut[0] || points_[i].getX() < cut[0]) {
					start_pos[0] = 0;
				}
				else if (points_[i - 1].getX() > cut[0] || points_[i].getX() > cut[0]) {
					start_pos[0] = 1;
				}
				else {
					start_pos[0] = -1;
				}
				if (points_[i - 1].getY() < cut[1] || points_[i].getY() < cut[1]) {
					start_pos[1] = 0;
				}
				else if (points_[i - 1].getY() > cut[1] || points_[i].getY() > cut[1]) {
					start_pos[1] = 1;
				}
				else {
					start_pos[1] = -1;
				}
			}
			else {
				sub_traj->points_.push_back(p);
				sub_traj->points_.push_back(points_[i]);
				if (p.getX() < cut[0] || points_[i].getX() < cut[0]) {
					start_pos[0] = 0;
				}
				else if (p.getX() > cut[0] || points_[i].getX() > cut[0]) {
					start_pos[0] = 1;
				}
				else {
					start_pos[0] = -1;
				}
				if (p.getY() < cut[1] || points_[i].getY() < cut[1]) {
					start_pos[1] = 0;
				}
				else if (p.getY() > cut[1] || points_[i].getY() > cut[1]) {
					start_pos[1] = 1;
				}
				else {
					start_pos[1] = -1;
				}
			}
			//cout << "start_pos: " << start_pos[0] << " " << start_pos[1] << endl;
		}
	}
	sub_area_id = 0;
	if (start_pos[0] == 1) {
		sub_area_id = 1;
	}
	if (start_pos[1] == 1) {
		sub_area_id += 2;
	}
	sub_trajectories[sub_area_id].push_back(sub_traj);
}

bool FYP::Trajectory::checkOverlap(Point& p1, Point& p2, double* query_max, double* query_min) {
	bool is_overlap1 = isLineRectangleOverlap(p1.getX(), p1.getY(), p2.getX(), p2.getY(), query_min[0], query_min[1], query_max[0], query_max[1]);
	bool is_overlap2 = isLineRectangleOverlap(p1.getX(), p1.getT(), p2.getX(), p2.getT(), query_min[0], query_min[2], query_max[0], query_max[2]);
	bool is_overlap3 = isLineRectangleOverlap(p1.getY(), p1.getT(), p2.getY(), p2.getT(), query_min[1], query_min[2], query_max[1], query_max[2]);
	// cout << "is_overlap: " << is_overlap1 << " " << is_overlap2 << " " << is_overlap3 << endl;
	bool is_overlap = is_overlap1 && is_overlap2 && is_overlap3;
	return is_overlap;
}

bool FYP::Trajectory::split(const unsigned dim, const double cut, std::vector<unsigned int>& split_pos, std::vector<Point>& split_point) {
	/*
	split a trajectory into a set of sub-trajctories based on the specified cut in the given dimension.
	*/
	double sp[3] = { 0.0, 0.0, 0.0 };
	bool has_split = false;
	for (unsigned int i = 0; i < points_.size() - 1; i++) {
		double max_value = std::max(points_[i].getValueFromDim(dim), points_[i + 1].getValueFromDim(dim));
		double min_value = std::min(points_[i].getValueFromDim(dim), points_[i + 1].getValueFromDim(dim));
		if (min_value <= cut && cut < max_value) {
			split_pos.push_back(i);
			double ratio = (cut - points_[i].getValueFromDim(dim)) / (points_[i + 1].getValueFromDim(dim) - points_[i].getValueFromDim(dim));
			for (unsigned int d = 0; d < 3; d++) {
				if (d == dim) {
					sp[d] = cut;
				}
				else {
					sp[d] = points_[i].getValueFromDim(d) + ratio * (points_[i + 1].getValueFromDim(d) - points_[i].getValueFromDim(d));
				}
			}
			split_point.emplace_back(sp);
			has_split = true;
		}
	}
	return has_split;

}

int FYP::Trajectory::compare(const unsigned int dim, const double cut) const {
	int points_below_cut = 0;
	int points_above_cut = 0;
	for (int i = 0; i < points_.size(); i++) {
		if (points_[i].getValueFromDim(dim) > cut) {
			points_above_cut += 1;
		}
		else if (points_[i].getValueFromDim(dim) < cut) {
			points_below_cut += 1;
		}
	}
	if (points_below_cut == 0) {
		return 1;
	}
	if (points_above_cut == 0) {
		return -1;
	}
	return 0;
}

int FYP::Trajectory::fastCompare(const unsigned int dim, const double cut) const {
	for (unsigned int i = 0; i < points_.size(); i++) {
		if (points_[i].getValueFromDim(dim) > cut) {
			return 1;
		}
		if (points_[i].getValueFromDim(dim) < cut) {
			return -1;
		}
	}
	return 0;
}

void FYP::Trajectory::clearPoints() {
	points_.clear();
}

bool FYP::Trajectory::isValid()const {
	return !points_.empty();
}

void FYP::Trajectory::print() const {
	std::cout << "trajectory " << id_ << ": [";
	for (const Point& p : points_) {
		std::cout << "(" << p.getX() << ", " << p.getY() << ", " << p.getT() << ") ";
	}
	std::cout << "]" << std::endl;
}

double FYP::Trajectory::getDistanceToPoint(Point& p) {
	double min_dist = DBL_MAX;
	for (int i = 0; i < points_.size() - 1; i++) {
		double d = getPDistance(p.getX(), p.getY(), points_[i].getX(), points_[i].getY(), points_[i + 1].getX(), points_[i + 1].getY());
		if (d < min_dist) {
			min_dist = d;
		}
	}
	return min_dist;
}

unsigned int FYP::Trajectory::serialize(char* content) {
	unsigned int pos = 0;
	memcpy(content + pos, &original_id_, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(content + pos, &id_, sizeof(id_));
	pos += sizeof(id_);
	unsigned int size = points_.size();
	memcpy(content + pos, &size, sizeof(size));
	pos += sizeof(size);
	for (int i = 0; i < points_.size(); i++) {
		memcpy(content + pos, points_[i].getDataAddr(), sizeof(double) * 3);
		pos += sizeof(double) * 3;
	}
	return pos;
}

unsigned int FYP::Trajectory::getByteSize() const {
	unsigned int size = (points_.size() - 1) * (sizeof(unsigned int) + sizeof(double) * 6);
	return size;
}

unsigned int FYP::Trajectory::deserialize(char* content) {
	unsigned int pos = 0;
	memcpy(&original_id_, content + pos, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(&id_, content + pos, sizeof(id_));
	pos += sizeof(id_);
	unsigned int size;
	memcpy(&size, content + pos, sizeof(size));
	pos += sizeof(size);
	points_ = std::vector<Point>(size);
	for (int i = 0; i < size; i++) {
		memcpy(points_[i].getDataAddr(), content + pos, sizeof(double) * 3);
		pos += sizeof(double) * 3;
	}
	return pos;
}

void FYP::Trajectory::fillSegment(std::vector<Segment*>&segments, int pos) {
	for (int i = 0; i < points_.size() - 1; i++) {
		segments[pos + i] = new Segment(original_id_, points_[i], points_[i + 1]);
	}
}
