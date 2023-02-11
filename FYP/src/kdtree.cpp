#include "kdtree.h"
//#define DEBUG
#undef DEBUG

#ifdef DEBUG
#define dlog(info) info
#else
#define dlog(info)
#endif

const double min_prec = 0.0001;

unsigned int GCOTrajNode_Segnum = 100;

double CrossProduct(double* v1, double* v2) {
	return v1[0] * v2[1] - v2[0] * v1[1];
}

bool LineRectangleOverlap(double x_start, double y_start, double x_end, double y_end, double left, double bottom, double right, double top) {
	//cout << x_start << " " << y_start << " " << x_end << " " << y_end << " " << top << " " << bottom << " " << left << " " << right << endl;
	if (left <= x_start && x_start <= right && bottom <= y_start && y_start <= top) {
		//the start point is inside the rectangle;
		return true;
	}
	if (left <= x_end && x_end <= right && bottom <= y_end && y_end <= top) {
		//the end point is inside the rectangle;
		return true;
	}
	double v1[2];
	double v2[2];
	double v3[2];
	double v4[2];
	double v5[2];
	double v6[2];
	v1[0] = x_start - right;
	v1[1] = y_start - bottom;
	v2[0] = left - right;
	v2[1] = top - bottom;
	v3[0] = x_end - right;
	v3[1] = y_end - bottom;

	v4[0] = left - x_end;
	v4[1] = top - y_end;
	v5[0] = x_start - x_end;
	v5[1] = y_start - y_end;
	v6[0] = right - x_end;
	v6[1] = bottom - y_end;

	/*cout << left << " " << right << " " << bottom << " " << top << endl;
	cout << x_start << " " << y_start << " " << x_end << " " << y_end << endl;
	cout << "cross product: " << CrossProduct(v1, v2) << " " << CrossProduct(v3, v2) << endl;
	*/if ((CrossProduct(v1, v2) * CrossProduct(v3, v2) <= 0) && (CrossProduct(v4, v5) * CrossProduct(v6, v5) <= 0)) {
	//segment intersect with the diagonal line1
		return true;
	}
	v1[0] = x_start - left;
	v1[1] = y_start - bottom;
	v2[0] = right - left;
	v2[1] = top - bottom;
	v3[0] = x_end - left;
	v3[1] = y_end - bottom;
	v4[0] = left - x_end;
	v4[1] = bottom - y_end;
	v5[0] = x_start - x_end;
	v5[1] = y_start - y_end;
	v6[0] = right - x_end;
	v6[1] = top - y_end;
	//cout << "cross product: " << CrossProduct(v1, v2) << " " << CrossProduct(v3, v2) << endl;
	if ((CrossProduct(v1, v2) * CrossProduct(v3, v2) <= 0) && (CrossProduct(v4, v5) * CrossProduct(v6, v5) <= 0)) {
		//segment intersect with the diagonal line2
		return true;
	}
	return false;
}

bool Kaiyu::Point::IsEqualTo(Point& point) const {
	if (point.values_[0] == values_[0] && point.values_[1] == values_[1] && point.values_[2] == values_[2]) {
		return true;
	}
	else {
		return false;
	}
}


void Kaiyu::Trajectory::AddPoint(const double x, const double y, const double t) {
	points_.emplace_back(x, y, t);
}
void Kaiyu::Trajectory::AddPoint(const Point& point) {
	points_.emplace_back(point.x(), point.y(), point.t());
}

bool Kaiyu::Trajectory::IsTemporalOverlap(double temporal_range[2]) {
	if (points_.front().t() < temporal_range[0] || points_.back().t() > temporal_range[1]) {
		return false;
	}
	return true;
}


bool Kaiyu::Trajectory::IsOverlap(const unsigned int dim, const double value) const {
	for (int i = 0; i < points_.size() - 1; i++) {
		double max_value = max(points_[i].value(dim), points_[i + 1].value(dim));
		double min_value = min(points_[i].value(dim), points_[i + 1].value(dim));
		if (min_value <= value && value < max_value) {
			return true;
		}
	}
	return false;
}



unsigned int Kaiyu::Trajectory::PageSize() const {
	unsigned int size = sizeof(double) * 3 * points_.size() + sizeof(unsigned int);
	return size;
}


int Kaiyu::Trajectory::FastCompare(const unsigned int dim, const double cut) const {
	for (unsigned int i = 0; i < points_.size(); i++) {
		if (points_[i].value(dim) > cut) {
			return 1;
		}
		if (points_[i].value(dim) < cut) {
			return -1;
		}
	}
	return 0;
}
int Kaiyu::Trajectory::Compare(const unsigned int dim, const double cut) const {
	int points_below_cut = 0;
	int points_above_cut = 0;
	for (int i = 0; i < points_.size(); i++) {
		if (points_[i].value(dim) > cut) {
			points_above_cut += 1;
		}
		else if (points_[i].value(dim) < cut) {
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

Kaiyu::Trajectory::Trajectory() {

}

Kaiyu::Trajectory::Trajectory(unsigned int id, unsigned int point_num) {
	original_id_ = id;
	id_ = id;
	points_ = vector<Point>(point_num);
}


bool Kaiyu::Trajectory::IsOverlap(double* query_max, double* query_min) {
	for (unsigned int i = 0; i < points_.size() - 1; i++) {
		bool is_overlap = CheckOverlap(points_[i], points_[i + 1], query_max, query_min);
		if (is_overlap) {
			return true;
		}
	}
	return false;
}

Kaiyu::Segment::Segment(double s[3], double t[3]) {
	start.Set(s[0], s[1], s[2]);
	end.Set(t[0], t[1], t[2]);
}

void Kaiyu::Segment::Set(Point& p1, Point& p2) {
	memcpy(start.data_addr(), p1.data_addr(), sizeof(double) * 3);
	memcpy(end.data_addr(), p2.data_addr(), sizeof(double) * 3);
}

bool Kaiyu::Segment::IsTemporalOverlap(double temporal_range[2]) {
	if (start.t() < temporal_range[0] || end.t() > temporal_range[1]) {
		return false;
	}
	return true;
}
bool Kaiyu::Segment::IsOverlap(double query_max[3], double query_min[3]) {
	bool is_overlap1 = LineRectangleOverlap(start.x(), start.y(), end.x(), end.y(), query_min[0], query_min[1], query_max[0], query_max[1]);
	bool is_overlap2 = LineRectangleOverlap(start.x(), start.t(), end.x(), end.t(), query_min[0], query_min[2], query_max[0], query_max[2]);
	bool is_overlap3 = LineRectangleOverlap(start.y(), start.t(), end.y(), end.t(), query_min[1], query_min[2], query_max[1], query_max[2]);
	bool is_overlap = is_overlap1 && is_overlap2 && is_overlap3;
	return is_overlap;
}

bool Kaiyu::Trajectory::CheckOverlap(Point& p1, Point& p2, double* query_max, double* query_min) {
	bool is_overlap1 = LineRectangleOverlap(p1.x(), p1.y(), p2.x(), p2.y(), query_min[0], query_min[1], query_max[0], query_max[1]);
	bool is_overlap2 = LineRectangleOverlap(p1.x(), p1.t(), p2.x(), p2.t(), query_min[0], query_min[2], query_max[0], query_max[2]);
	bool is_overlap3 = LineRectangleOverlap(p1.y(), p1.t(), p2.y(), p2.t(), query_min[1], query_min[2], query_max[1], query_max[2]);
	// cout << "is_overlap: " << is_overlap1 << " " << is_overlap2 << " " << is_overlap3 << endl;
	bool is_overlap = is_overlap1 && is_overlap2 && is_overlap3;
	return is_overlap;
}

void Kaiyu::Trajectory::Clear() {
	points_.clear();
}

void Kaiyu::Trajectory::Print() const {
	cout << "trajectory " << id_ << ": [";
	for (const Point& p : points_) {
		cout << "(" << p.x() << ", " << p.y() << ", " << p.t() << ") ";
	}
	cout << "]" << endl;
}

bool Kaiyu::Trajectory::IsValid()const {
	return !points_.empty();
}

void Kaiyu::Trajectory::CrossSplit(const double cut[2], vector<list<Trajectory*> >& sub_trajectories) {
	//cout << "Cross splitting trajectory with" << cut[0]<<" "<<cut[1]<<endl;
	//Print();
	int start_pos[2] = { -1, -1 };
	if (points_[0].x() < cut[0]) {
		start_pos[0] = 0;
	}
	if (points_[0].x() > cut[0]) {
		start_pos[0] = 1;
	}
	if (points_[0].y() < cut[1]) {
		start_pos[1] = 0;
	}
	if (points_[0].y() > cut[1]) {
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
		if (points_[i].x() < cut[0]) {
			next_pos[0] = 0;
		}
		else if (points_[i].x() > cut[0]) {
			next_pos[0] = 1;
		}
		else {
			next_pos[0] = -1;
		}

		if (points_[i].y() < cut[1]) {
			next_pos[1] = 0;
		}
		else if (points_[i].y() > cut[1]) {
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
			double ratio1 = (cut[0] - points_[i - 1].x()) / (points_[i].x() - points_[i - 1].x());
			double ratio2 = (cut[1] - points_[i - 1].y()) / (points_[i].y() - points_[i - 1].y());
			double small_ratio = min(ratio1, ratio2);
			double large_ratio = max(ratio1, ratio2);
			//cout << "ratio: " << small_ratio << " " << large_ratio << endl;
			//cout << (small_ratio == 0) << " " << (large_ratio == 0) << endl;
			double coord1[3];
			double coord2[3];
			coord1[0] = points_[i - 1].x() + small_ratio * (points_[i].x() - points_[i - 1].x());
			coord1[1] = points_[i - 1].y() + small_ratio * (points_[i].y() - points_[i - 1].y());
			coord1[2] = points_[i - 1].t() + small_ratio * (points_[i].t() - points_[i - 1].t());
			coord2[0] = points_[i - 1].x() + large_ratio * (points_[i].x() - points_[i - 1].x());
			coord2[1] = points_[i - 1].y() + large_ratio * (points_[i].y() - points_[i - 1].y());
			coord2[2] = points_[i - 1].t() + large_ratio * (points_[i].t() - points_[i - 1].t());
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
					if (0.5 * (p1.x() + p2.x()) > cut[0]) {
						sub_area_id = 1;
					}
					if (0.5 * (p1.y() + p2.y()) > cut[1]) {
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
				if (0.5 * (p1.x() + p2.x()) > cut[0]) {
					sub_area_id = 1;
				}
				if (0.5 * (p1.y() + p2.y()) > cut[1]) {
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
			double ratio = (cut[dim] - points_[i - 1].value(dim)) / (points_[i].value(dim) - points_[i - 1].value(dim));
			//cout << "ratio: " << ratio << endl;
			double coord[3];
			coord[0] = points_[i - 1].x() + ratio * (points_[i].x() - points_[i - 1].x());
			coord[1] = points_[i - 1].y() + ratio * (points_[i].y() - points_[i - 1].y());
			coord[2] = points_[i - 1].t() + ratio * (points_[i].t() - points_[i - 1].t());
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
				if (points_[i - 1].x() < cut[0] || points_[i].x() < cut[0]) {
					start_pos[0] = 0;
				}
				else if (points_[i - 1].x() > cut[0] || points_[i].x() > cut[0]) {
					start_pos[0] = 1;
				}
				else {
					start_pos[0] = -1;
				}
				if (points_[i - 1].y() < cut[1] || points_[i].y() < cut[1]) {
					start_pos[1] = 0;
				}
				else if (points_[i - 1].y() > cut[1] || points_[i].y() > cut[1]) {
					start_pos[1] = 1;
				}
				else {
					start_pos[1] = -1;
				}
			}
			else {
				sub_traj->points_.push_back(p);
				sub_traj->points_.push_back(points_[i]);
				if (p.x() < cut[0] || points_[i].x() < cut[0]) {
					start_pos[0] = 0;
				}
				else if (p.x() > cut[0] || points_[i].x() > cut[0]) {
					start_pos[0] = 1;
				}
				else {
					start_pos[0] = -1;
				}
				if (p.y() < cut[1] || points_[i].y() < cut[1]) {
					start_pos[1] = 0;
				}
				else if (p.y() > cut[1] || points_[i].y() > cut[1]) {
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

bool Kaiyu::Trajectory::Split(const unsigned dim, const double cut, vector<unsigned int>& split_pos, vector<Point>& split_point) {
	/*
	split a trajectory into a set of sub-trajctories based on the specified cut in the given dimension.
	*/
	double sp[3] = { 0.0, 0.0, 0.0 };
	bool has_split = false;
	for (unsigned int i = 0; i < points_.size() - 1; i++) {
		double max_value = max(points_[i].value(dim), points_[i + 1].value(dim));
		double min_value = min(points_[i].value(dim), points_[i + 1].value(dim));
		if (min_value <= cut && cut < max_value) {
			split_pos.push_back(i);
			double ratio = (cut - points_[i].value(dim)) / (points_[i + 1].value(dim) - points_[i].value(dim));
			for (unsigned int d = 0; d < 3; d++) {
				if (d == dim) {
					sp[d] = cut;
				}
				else {
					sp[d] = points_[i].value(d) + ratio * (points_[i + 1].value(d) - points_[i].value(d));
				}
			}
			split_point.emplace_back(sp);
			has_split = true;
		}
	}
	return has_split;

}

double Kaiyu::GCOTrajNode::L2Distance(Point& p) {
	if (min_range[0] < p.x() && p.x() < max_range[0]) {
		if (p.y() > max_range[1]) {
			return p.y() - max_range[1];
		}
		if (p.y() < min_range[1]) {
			return min_range[1] - p.y();
		}
		return 0;
	}
	else {
		if (min_range[1] < p.y() && p.y() < max_range[1]) {
			if (p.x() > max_range[0]) {
				return p.x() - max_range[0];
			}
			else {
				return min_range[0] - p.x();
			}
		}
		double dist = 0;
		double anchor[2];
		if (p.x() > max_range[0]) {
			anchor[0] = max_range[0];
		}
		else {
			anchor[0] = min_range[0];
		}
		if (p.y() > max_range[1]) {
			anchor[1] = max_range[1];
		}
		else {
			anchor[1] = min_range[1];
		}
		dist = sqrt((anchor[0] - p.x()) * (anchor[0] - p.x()) + (anchor[1] - p.y()) * (anchor[1] - p.y()));
		return dist;
	}

}

double Kaiyu::TrajStoreNode::L2Distance(Point& p) {
	if (min_range[0] < p.x() && p.x() < max_range[0]) {
		if (p.y() > max_range[1]) {
			return p.y() - max_range[1];
		}
		if (p.y() < min_range[1]) {
			return min_range[1] - p.y();
		}
		return 0;
	}
	else {
		if (min_range[1] < p.y() && p.y() < max_range[1]) {
			if (p.x() > max_range[0]) {
				return p.x() - max_range[0];
			}
			else {
				return min_range[0] - p.x();
			}
		}
		double dist = 0;
		double anchor[2];
		if (p.x() > max_range[0]) {
			anchor[0] = max_range[0];
		}
		else {
			anchor[0] = min_range[0];
		}
		if (p.y() > max_range[1]) {
			anchor[1] = max_range[1];
		}
		else {
			anchor[1] = min_range[1];
		}
		dist = sqrt((anchor[0] - p.x()) * (anchor[0] - p.x()) + (anchor[1] - p.y()) * (anchor[1] - p.y()));
		return dist;
	}
}


bool Kaiyu::KDTreeNode::IsTemporalOverlap(double temporal_range[2]) {
	if (max_range[2] < temporal_range[0] || min_range[2] > temporal_range[1]) {
		return false;
	}
	return true;
}

bool Kaiyu::TrajStoreNode::IsTemporalOverlap(const double temporal_interval[2]) {
	if (max_range[2] < temporal_interval[0] || min_range[2] > temporal_interval[1]) {
		return false;
	}
	return true;
}

bool Kaiyu::TrajStoreNode::IsOverlap(const double query_max[3], const double query_min[3]) {
	for (int i = 0; i < 3; i++) {
		if (max_range[i] < query_min[i] || min_range[i] > query_max[i]) {
			return false;
		}
	}
	return true;
}

bool Kaiyu::GCOTrajNode::IsOverlap(const double query_max[3], const double query_min[3]) {
	for (int i = 0; i < 3; i++) {
		if (max_range[i] < query_min[i] || min_range[i] > query_max[i]) {
			return false;
		}
	}
	return true;
}

bool Kaiyu::TBTreeNode::IsOverlap(const double range_max[3], const double range_min[3]) {
	for (int i = 0; i < 3; i++) {
		if (range_max[i] < min_range[i] || range_min[i] > max_range[i]) {
			return false;
		}
	}
	return true;
}

bool Kaiyu::TBTreeNode::IsTemporalOverlap(double temporal_range[2]) {
	if (min_range[2] > temporal_range[1] || max_range[2] < temporal_range[0]) {
		return false;
	}
	return true;
}

bool Kaiyu::GCOTrajNode::IsTemporalOverlap(const double temporal_range[2]) {
	if (min_range[2] > temporal_range[1] || max_range[2] < temporal_range[0]) {
		return false;
	}
	return true;
}

bool Kaiyu::KDTreeNode::IsOverlap(double* query_max, double* query_min) {
	for (int i = 0; i < 3; i++) {
		if (query_max[i] < min_range[i] || query_min[i] > max_range[i]) {
			return false;
		}
	}
	return true;
}

void Kaiyu::KDTree::ComputeSegment(const unsigned int id) {
	KDTreeNode* node = tree_nodes_[id];
	node->segment_num_ = 0;
	for (int id : node->data_) {
		node->segment_num_ += trajectories_[id]->SegmentNum();
	}
}

void Kaiyu::KDTreeNode::Print() const {
	cout << "[";
	for (int i = 0; i < 3; i++) {
		cout << "(" << min_range[i] << ", " << max_range[i] << ") ";
	}
	cout << "]" << endl;
}

void Kaiyu::KDTree::ComputeMBRange(const unsigned int id) {
	KDTreeNode* node = tree_nodes_[id];
	node->idx_range_max[0] = 0;
	node->idx_range_max[1] = 0;
	node->idx_range_max[2] = 0;
	node->idx_range_max[0] = cuts_[0].size();
	node->idx_range_max[1] = cuts_[1].size();
	node->idx_range_max[2] = cuts_[2].size();
	for (int dim = 0; dim < 3; dim++) {
		for (int i = 0; i < cuts_[dim].size(); i++) {
			if (cuts_[dim][i] >= node->max_range[dim]) {
				node->idx_range_max[dim] = i;
				break;
			}
		}
		for (int i = cuts_[dim].size() - 1; i >= 0; i--) {
			if (cuts_[dim][i] <= node->min_range[dim]) {
				node->idx_range_min[dim] = i;
				break;
			}
		}
	}
}

void Kaiyu::KDTree::ComputeMBR(const unsigned int id) {
	KDTreeNode* node = tree_nodes_[id];
	node->max_range[0] = -DBL_MAX;
	node->max_range[1] = -DBL_MAX;
	node->max_range[2] = -DBL_MAX;
	node->min_range[0] = DBL_MAX;
	node->min_range[1] = DBL_MAX;
	node->min_range[2] = DBL_MAX;
	for (int id : node->data_) {
		Trajectory* traj = trajectories_[id];
		for (auto const& p : traj->points_) {
			node->max_range[0] = max(node->max_range[0], p.x());
			node->max_range[1] = max(node->max_range[1], p.y());
			node->max_range[2] = max(node->max_range[2], p.t());
			node->min_range[0] = min(node->min_range[0], p.x());
			node->min_range[1] = min(node->min_range[1], p.y());
			node->min_range[2] = min(node->min_range[2], p.t());
		}
	}
}

bool Kaiyu::KDTree::IsSplitable(const unsigned int node_id, const unsigned int dim) {
	KDTreeNode* node = tree_nodes_[node_id];
	if (node->max_range[dim] - node->min_range[dim] < min_prec) {
		return false;
	}
	else {
		return true;
	}
	//for (int traj_id : node->data_) {
	//	Trajectory* traj = trajectories_[traj_id];
	//	for (Point& p : traj->points_) {
	//		values.insert(p.value(dim));
	//		if (values.size() > 1) {
	//			return true;
	//		}
	//	}
	//}
	//return false;
}

void Kaiyu::KDTree::SplitTrajectory(const unsigned int traj_id, const unsigned int dim, const double cut, list<int>& sub_traj1, list<int>& sub_traj2) {
	Trajectory* traj = trajectories_[traj_id];
	list<Trajectory*> subtrajectories;
	subtrajectories.push_back(new Trajectory());
	Trajectory* ptr = subtrajectories.back();
	ptr->original_id_ = traj->original_id_;
	ptr->id_ = traj->id_;
	ptr->AddPoint(traj->points_[0]);
	int is_in_side = 0;
	if (traj->points_[0].value(dim) < cut) {
		is_in_side = 1;
	}
	else if (traj->points_[0].value(dim) > cut) {
		is_in_side = 2;
	}
	double mid_point[3];
	int num = 0;
	for (int i = 1; i < traj->points_.size(); i++) {
		int side = 0;
		if (traj->points_[i].value(dim) < cut) {
			side = 1;
		}
		else if (traj->points_[i].value(dim) > cut) {
			side = 2;
		}
		if (side + is_in_side == 3) {
			if (is_in_side == 1) {
				sub_traj1.push_back(ptr->id_);
			}
			else {
				sub_traj2.push_back(ptr->id_);
			}
			if (traj->points_[i - 1].value(dim) == cut) {
				subtrajectories.push_back(new Trajectory());
				ptr = subtrajectories.back();
				ptr->AddPoint(traj->points_[i - 1]);
				ptr->original_id_ = traj->original_id_;
				ptr->id_ = trajectories_.size() + num;
				num += 1;
			}
			else {
				double ratio = (cut - traj->points_[i - 1].value(dim)) / (traj->points_[i].value(dim) - traj->points_[i - 1].value(dim));
				for (int d = 0; d < 3; d++) {
					mid_point[d] = traj->points_[i - 1].value(d) + ratio * (traj->points_[i].value(d) - traj->points_[i - 1].value(d));
				}
				ptr->AddPoint(mid_point);
				subtrajectories.push_back(new Trajectory());
				ptr = subtrajectories.back();
				ptr->original_id_ = traj->original_id_;
				ptr->id_ = trajectories_.size() + num;
				num += 1;
				ptr->AddPoint(mid_point);
			}
			is_in_side = side;
		}
		else {
			if (side > 0 && is_in_side == 0) {
				is_in_side = side;
			}
		}
		ptr->AddPoint(traj->points_[i]);

	}
	if (is_in_side == 1) {
		sub_traj1.push_back(ptr->id_);
	}
	else {
		sub_traj2.push_back(ptr->id_);
	}
	for (Trajectory* ptr : subtrajectories) {
		if (traj_id == ptr->id_) {
			delete trajectories_[traj_id];
			trajectories_[traj_id] = ptr;
		}
		else {
			trajectories_.push_back(ptr);
		}
	}
	subtrajectories.clear();
}
/*
void Kaiyu::KDTree::SplitTrajectory(const unsigned int traj_id, const unsigned int dim, const double cut, list<int>& sub_traj1, list<int>& sub_traj2) {
	Trajectory* traj = trajectories_[traj_id];
	vector<unsigned int> locs;
	vector<Point> points;
	traj->Split(dim, cut, locs, points);
	if (locs.size() > 0) {
		list<Trajectory*> sub_trajectories;
		unsigned int point_id = 0;
		unsigned int split_loc_id = 0;
		sub_trajectories.push_back(nullptr);
		sub_trajectories.back() = new Trajectory();
		sub_trajectories.back()->original_id_ = traj->original_id_;
		while (point_id < traj->points_.size()) {
			if (split_loc_id == locs.size() || point_id < locs[split_loc_id]) {
				sub_trajectories.back()->AddPoint(traj->points_[point_id]);
			}
			else if (point_id == locs[split_loc_id]) {
				sub_trajectories.back()->AddPoint(traj->points_[point_id]);
				bool cut_in_the_middle = !points[split_loc_id].IsEqualTo(traj->points_[point_id]);
				if(cut_in_the_middle){
					sub_trajectories.back()->AddPoint(points[split_loc_id]);
				}
				sub_trajectories.push_back(nullptr);
				sub_trajectories.back() = new Trajectory();
				sub_trajectories.back()->original_id_ = traj->original_id_;
				if (cut_in_the_middle) {
					sub_trajectories.back()->AddPoint(points[split_loc_id]);
				}
				else {
					sub_trajectories.back()->AddPoint(traj->points_[point_id]);
				}
				split_loc_id += 1;

			}
			point_id += 1;
		}
		auto traj_ptr = sub_trajectories.begin();
		(*traj_ptr)->id_ = traj_id;
		delete trajectories_[traj_id];
		trajectories_[traj_id] = *traj_ptr;
		++traj_ptr;
		while (traj_ptr != sub_trajectories.end()) {
			(*traj_ptr)->id_ = trajectories_.size();
			trajectories_.push_back(*traj_ptr);
			++traj_ptr;
		}
		for (auto traj_ptr : sub_trajectories) {
			if (traj_ptr->FastCompare(dim, cut) <= 0) {
				sub_traj1.push_back(traj_ptr->id_);
			}
			else {
				sub_traj2.push_back(traj_ptr->id_);
			}
		}
	}
}

*/
void Kaiyu::KDTree::Initialize(const string traj_file_name, const string cut_file) {
	std::ifstream ifs(traj_file_name, std::ios::in);
	unsigned int total_traj = 0;
	unsigned int point_num = 0;
	ifs >> total_traj;
	original_trajectories_ = vector<Trajectory*>(total_traj, nullptr);
	trajectories_ = vector<Trajectory*>(total_traj, nullptr);
	for (unsigned int i = 0; i < total_traj; i++) {
		ifs >> point_num;
		original_trajectories_[i] = new Trajectory(i, point_num);
		trajectories_[i] = new Trajectory(i, point_num);;
		for (unsigned int point_id = 0; point_id < point_num; point_id += 1) {
			double x, y, t;
			ifs >> x >> y >> t;
			original_trajectories_[i]->points_[point_id].Set(x, y, t);
			trajectories_[i]->points_[point_id].Set(x, y, t);
		}
	}
	ifs.close();
	cout << "trajectory size: " << trajectories_.size() << endl;
	ifs.open(cut_file);
	int cut_num = 0;
	ifs >> cut_num;
	cuts_ = vector<vector<double> >(3);
	for (unsigned int i = 0; i < 3; i++) {
		cuts_[i] = vector<double>(cut_num);
		for (unsigned int cut_id = 0; cut_id < cut_num; cut_id++) {
			ifs >> cuts_[i][cut_id];
		}
	}
	ifs.close();
	tree_nodes_.push_back(nullptr);
	tree_nodes_[0] = new KDTreeNode();
	tree_nodes_[0]->depth_ = 0;
	tree_nodes_[0]->data_ = vector<int>(total_traj);
	for (unsigned int i = 0; i < total_traj; i++) {
		tree_nodes_[0]->data_[i] = i;
	}
	tree_nodes_[0]->children[0] = 0;
	tree_nodes_[0]->children[1] = 0;
	ComputeMBR(0);
	ComputeMBRange(0);
	ComputeSegment(0);
	tree_nodes_[0]->overall_segment_ = 0;
	tree_nodes_[0]->is_leaf_ = true;
	depth_ = 0;



	initialized_ = true;
}

Kaiyu::KDTree::KDTree() {
	propagated_ = false;
	initialized_ = false;
}

Kaiyu::TrajStore::TrajStore() {
	page_size_ = 1024;
	stop_increase_ratio = 0;
	node_access = 0;
	leaf_access = 0;
	block_acces = 0;
}

Kaiyu::TrajStore::~TrajStore() {
	for (TrajStoreNode* node : tree_nodes_) {
		delete node;
	}
	tree_nodes_.clear();
}

Kaiyu::TrajStoreNode::TrajStoreNode() {
	for (int i = 0; i < 3; i++) {
		max_range[i] = -DBL_MAX;
		min_range[i] = DBL_MAX;
	}
	children_ = nullptr;
	data_ = nullptr;
	content = nullptr;
	id_ = -1;
	is_leaf_ = true;
}

void Kaiyu::TrajStoreNode::UpdateRange(double x, double y, double t) {
	max_range[0] = max(max_range[0], x);
	min_range[0] = min(min_range[0], x);
	max_range[1] = max(max_range[1], y);
	min_range[1] = min(min_range[1], y);
	max_range[2] = max(max_range[2], t);
	min_range[2] = min(min_range[2], t);
}

void Kaiyu::TrajStore::Initialize(const string traj_file_name) {
	std::ifstream ifs(traj_file_name, std::ios::in);
	unsigned int total_traj = 0;
	unsigned int point_num = 0;
	ifs >> total_traj;
	cout << total_traj << " trajectories" << endl;
	TrajStoreNode* root = new TrajStoreNode();
	root->data_ = new vector<Trajectory*>(total_traj, nullptr);
	root->id_ = tree_nodes_.size();
	cout << "root node constructed" << endl;
	for (unsigned int i = 0; i < total_traj; i++) {
		ifs >> point_num;
		//cout << i << " " << point_num << endl;
		root->data_->at(i) = new Trajectory(i, point_num);
		for (unsigned int point_id = 0; point_id < point_num; point_id += 1) {
			double x, y, t;
			ifs >> x >> y >> t;
			root->UpdateRange(x, y, t);
			root->data_->at(i)->points_[point_id].Set(x, y, t);
		}
	}
	ifs.close();
	cout << "trajectories loaded" << endl;
	getchar();
	tree_nodes_.push_back(root);
	list<TrajStoreNode*> queue;
	queue.push_back(root);
	while (!queue.empty()) {
		TrajStoreNode* node = queue.front();
		queue.pop_front();
		bool success = StopCondition(node->id_);

		if (success) {
			for (int child_id : *(node->children_)) {
				TrajStoreNode* child = tree_nodes_[child_id];
				queue.push_back(child);
			}
		}
		else {
			FinalizeNode(node->id_);
		}
		//cout <<"success: "<< success << " queue.size(): " << queue.size() << endl;
	}
}

void Kaiyu::Trajectory::FillSegment(vector<Segment*>& segments, int pos) {
	for (int i = 0; i < points_.size() - 1; i++) {
		segments[pos + i] = new Segment();
		segments[pos + i]->original_id_ = original_id_;
		memcpy(segments[pos + i]->start.data_addr(), points_[i].data_addr(), sizeof(double) * 3);
		memcpy(segments[pos + i]->end.data_addr(), points_[i + 1].data_addr(), sizeof(double) * 3);
	}
}

unsigned int Kaiyu::Segment::Serialization(char* content) {
	unsigned int pos = 0;
	memcpy(content + pos, &original_id_, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(content + pos, start.data_addr(), sizeof(double) * 3);
	pos += sizeof(double) * 3;
	memcpy(content + pos, end.data_addr(), sizeof(double) * 3);
	pos += sizeof(double) * 3;
	return pos;
}

unsigned int Kaiyu::Segment::Deserialization(char* content) {
	unsigned int pos = 0;
	memcpy(&original_id_, content + pos, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(start.data_addr(), content + pos, sizeof(double) * 3);
	pos += sizeof(double) * 3;
	memcpy(end.data_addr(), content + pos, sizeof(double) * 3);
	pos += sizeof(double) * 3;
	return pos;
}

unsigned int Kaiyu::Segment::ByteSize() const {
	unsigned int size = 0;
	size += sizeof(original_id_);
	size += sizeof(double) * 6;
	return size;
}

void Kaiyu::TrajStore::FinalizeNode(int node_id) {
	TrajStoreNode* node = tree_nodes_[node_id];
	if (node->data_ == nullptr) {
		cout << "Node " << node_id << " does not have any data yet." << endl;
		exit(0);
	}
	int total_segments_num = 0;
	for (Trajectory* traj : *(node->data_)) {
		total_segments_num += traj->points_.size() - 1;
	}
	node->is_leaf_ = false;
	vector<Segment*> all_segments(total_segments_num, nullptr);
	int start_pos = 0;
	for (Trajectory* traj : *(node->data_)) {
		traj->FillSegment(all_segments, start_pos);
		start_pos += traj->points_.size() - 1;
	}
	sort(all_segments.begin(), all_segments.end(), [](Segment* lhs, Segment* rhs) {return lhs->end.t() < rhs->end.t(); });
	unsigned int size = all_segments.size() * (sizeof(unsigned int) + sizeof(double) * 3 * 2);
	int child_num = (int)ceil(1.0 * all_segments.size() / B);
	node->children_ = new vector<int>(child_num);
	for (int c_id = 0; c_id < child_num; c_id++) {

		TrajStoreNode* child = new TrajStoreNode();
		child->id_ = tree_nodes_.size();
		child->is_leaf_ = true;
		node->children_->at(c_id) = child->id_;
		child->seg_num_ = all_segments.size() - c_id * B;
		if (child->seg_num_ > B) {
			child->seg_num_ = B;
		}
		child->content = new char[child->seg_num_ * (sizeof(unsigned int) + sizeof(double) * 6)];
		unsigned int start_idx = 0;
		for (int i = 0; i < child->seg_num_; i++) {
			Segment* seg = all_segments[c_id * B + i];
			child->UpdateRange(seg->start.x(), seg->start.y(), seg->start.t());
			child->UpdateRange(seg->end.x(), seg->end.y(), seg->end.t());
			int len = seg->Serialization(child->content + start_idx);
			start_idx += len;
		}
		tree_nodes_.push_back(child);
	}
	for (Trajectory* traj : *(node->data_)) {
		delete traj;
	}
	delete node->data_;
	node->data_ = nullptr;
}



bool Kaiyu::TrajStore::StopCondition(int node_id) {
	TrajStoreNode* node = tree_nodes_[node_id];
	//cout << "node range: " << node->min_range[0] << " " << node->max_range[0] << ", " << node->min_range[1] << " " << node->max_range[1] << endl;
	unsigned int byte_size = 0;
	unsigned int seg_num_before_split = 0;
	unsigned int seg_num_after_split = 0;
	//cout << "trajectory: " << endl;
	for (int i = 0; i < node->data_->size(); i++) {
		byte_size += node->data_->at(i)->ByteSize();
		//node->data_->at(i)->Print();
		seg_num_before_split += node->data_->at(i)->SegmentNum();
	}
	if (byte_size <= page_size_) {
		//no need to split
		return false;
	}
	double cost_not_split = (query_spatial_domain[0] + node->x_range()) * (query_spatial_domain[1] + node->y_range()) * ceil(1.0 * byte_size / page_size_);

	TrajStoreNode* children[4];
	for (int i = 0; i < 4; i++) {
		children[i] = new TrajStoreNode();
		children[i]->data_ = new vector<Trajectory*>();
		children[i]->data_->reserve(node->data_->size());
	}
	vector<list<Trajectory*> > sub_trajectories(4);
	double cross_cuts[2] = { 0.5 * (node->max_range[0] + node->min_range[0]), 0.5 * (node->max_range[1] + node->min_range[1]) };
	for (int traj_id = 0; traj_id < node->data_->size(); traj_id++) {
		Trajectory* traj = node->data_->at(traj_id);
		traj->CrossSplit(cross_cuts, sub_trajectories);
		for (int child_id = 0; child_id < 4; child_id++) {
			for (Trajectory* subtraj : sub_trajectories[child_id]) {
				//subtraj->Print();
				children[child_id]->data_->push_back(subtraj);
				seg_num_after_split += subtraj->SegmentNum();
			}
			sub_trajectories[child_id].clear();
		}
	}
	double cost_after_split = 0;
	unsigned int split_byte_sizes[4] = { 0, 0, 0, 0 };
	for (int child_id = 0; child_id < 4; child_id++) {
		for (Trajectory* traj : *(children[child_id]->data_)) {
			split_byte_sizes[child_id] += traj->ByteSize();
			for (Point& p : traj->points_) {
				children[child_id]->UpdateRange(p.x(), p.y(), p.t());
			}
		}
		cost_after_split += (query_spatial_domain[0] + 0.5 * node->x_range()) * (query_spatial_domain[1] + 0.5 * node->y_range()) * ceil(1.0 * split_byte_sizes[child_id] / page_size_);
	}
	/*cout << "query_spatial domain: " << query_spatial_domain[0] << " " << query_spatial_domain[1] << endl;
	cout << "x_range and y_range: " << node->x_range() << " " << node->y_range() << endl;
	cout << "p/page_size: " << ceil(1.0 * byte_size / page_size_) << endl;
	cout << sizeof(unsigned int) + sizeof(double) * 6 << endl;
	cout << byte_size << " ";
	for (int i = 0; i < 4; i++) {
		cout << split_byte_sizes[i] << " ";
	}
	cout << endl;
	cout << "cost not split: " << cost_not_split << ", cost_after_split: " << cost_after_split << endl;
	cout << "seg num before: " << seg_num_before_split << ", after: " << seg_num_after_split << endl;
	cout << "four children: " << endl;*/
	/*for (int i = 0; i < 4; i++) {
		cout << "children " << i << ": " << children[i]->min_range[0] << " " << children[i]->max_range[0] << ", " << children[i]->min_range[1] << " " << children[i]->max_range[1] << endl;
		for (int i = 0; i < children[i]->data_->size(); i++) {
			children[i]->data_->at(i)->Print();
		}
		cout << "####" << endl;
	}*/
	if (seg_num_before_split * (1 + stop_increase_ratio) < seg_num_after_split) {
		for (int child_id = 0; child_id < 4; child_id++) {
			for (Trajectory* traj : *(children[child_id]->data_)) {
				delete traj;
			}
			delete children[child_id]->data_;
			delete children[child_id];
		}
		return false;
	}
	if (cost_not_split > (1 + epsilon) * cost_after_split) {
		node->children_ = new vector<int>();
		for (int i = 0; i < 4; i++) {
			if (children[i]->data_->empty()) {

				if (children[i]->data_ == nullptr) {
					cout << "null" << endl;
				}
				if (children[i] == nullptr) {
					cout << "null2" << endl;
				}
				delete children[i]->data_;
				delete children[i];
			}
			else {
				children[i]->id_ = tree_nodes_.size();
				tree_nodes_.push_back(children[i]);
				node->children_->push_back(children[i]->id_);
			}

		}
		for (int i = 0; i < node->data_->size(); i++) {
			delete node->data_->at(i);
		}
		delete node->data_;
		node->data_ = nullptr;
		return true;
	}
	else {
		for (int child_id = 0; child_id < 4; child_id++) {
			for (Trajectory* traj : *(children[child_id]->data_)) {
				delete traj;
			}
			delete children[child_id]->data_;
			delete children[child_id];
		}
		return false;
	}
}


void Kaiyu::KDTree::KNNQuery(int k, double temporal_range[2], double coord[2]) {
	node_access_ = 0;
	leaf_node_access_ = 0;
	KDTreeNode* root = tree_nodes_[0];
	Point q;
	q.Set(coord[0], coord[1], 0);
	priority_queue<pair<double, int>, vector<pair<double, int> >, std::greater<pair<double, int> > > pqueue;
	if (root->IsTemporalOverlap(temporal_range)) {
		double dist = root->L2Distance(q);
		pqueue.emplace(dist, 0);
	}
	vector<pair<double, int> > results;
	while (!pqueue.empty()) {
		if (results.size() == k && results.back().first < pqueue.top().first) {
			break;
		}
		pair<double, int> top = pqueue.top();
		pqueue.pop();
		KDTreeNode* node = tree_nodes_[top.second];
		if (node->is_leaf_) {
			leaf_node_access_ += 1;
			node_access_ += 1;
			for (int traj_id : node->data_) {
				Trajectory* traj = trajectories_[traj_id];
				if (!traj->IsTemporalOverlap(temporal_range)) {
					continue;
				}
				double distance = traj->DistanceToPoint(q);
				int original_id = traj->original_id_;
				bool is_in = false;
				bool changed = false;
				for (int i = 0; i < results.size(); i++) {
					if (results[i].second == original_id) {
						is_in = true;
						if (results[i].first > distance) {
							results[i].first = distance;
							changed = true;
						}
						break;
					}
				}
				if (!is_in) {
					if (results.size() < k) {
						results.emplace_back(distance, original_id);
						changed = true;
					}
					else {
						if (results.back().first > distance) {
							results.back().first = distance;
							results.back().second = original_id;
							changed = true;
						}
					}
				}
				if (changed) {
					sort(results.begin(), results.end());
				}
			}
		}
		else {
			node_access_ += 1;
			for (int i = 0; i < 2; i++) {
				int child_id = node->children[i];
				KDTreeNode* child = tree_nodes_[child_id];
				if (child->IsTemporalOverlap(temporal_range)) {
					double dist = child->L2Distance(q);
					pqueue.emplace(dist, child_id);
				}
			}
		}
	}
}

Kaiyu::GCOTrajNode::GCOTrajNode(int id) {
	id_ = id;
	seg_num_ = 0;
	for (int i = 0; i < 3; i++) {
		max_range[i] = -DBL_MAX;
		min_range[i] = DBL_MAX;
	}
}
Kaiyu::GCOTrajNode::GCOTrajNode(double max_r[3], double min_r[3], int id) {
	for (int i = 0; i < 3; i++) {
		max_range[i] = max_r[i];
		min_range[i] = min_r[i];
	}
	seg_num_ = 0;
	id_ = id;
}

void Kaiyu::GCOTrajNode::UpdateRange(double value[3]) {
	for (int i = 0; i < 3; i++) {
		max_range[i] = max(max_range[i], value[i]);
		min_range[i] = min(min_range[i], value[i]);
	}
}

void Kaiyu::GCOTrajNode::AddSegment(Segment& seg) {
	seg_num_ += 1;
	cache.emplace_back(seg.start.data_addr(), seg.end.data_addr());
	cache.back().original_id_ = seg.original_id_;
	if (cache.size() == GCOTrajNode_Segnum) {
		WritePage();
		cache.clear();
	}
}

void Kaiyu::GCOTrajNode::WritePage() {
	int pos = pages_.size();
	unsigned int seg_num = cache.size(); 
	pages_.push_back(nullptr);
	unsigned int page_size = sizeof(unsigned int) + cache[0].ByteSize() * seg_num;
	pages_[pos] = new char[page_size];
	unsigned int page_pos = 0;
	memcpy(pages_[pos] + page_pos, &seg_num, sizeof(unsigned int));
	page_pos += sizeof(unsigned int);
	for (unsigned int i = 0; i < seg_num; i++) {
		cache[i].Serialization(pages_[pos] + page_pos);
		page_pos += cache[i].ByteSize();
	}
	cache.clear();
}

void Kaiyu::GCOTrajNode::LoadPage(vector<Segment>& segments, unsigned int page_id) {
	unsigned int seg_num;
	memcpy(&seg_num, pages_[page_id], sizeof(unsigned int));
	segments.resize(seg_num);
	unsigned int pos = sizeof(unsigned int);
	for (unsigned int i = 0; i < seg_num; i++) {
		segments[i].Deserialization(pages_[page_id] + pos);
		pos += segments[i].ByteSize();
	}
}

Kaiyu::GCOTraj::GCOTraj() {
	leaf_access = 0;
	node_access = 0;
	block_access = 0;
	grid_granularity[0] = grid_granularity[1] = grid_granularity[2] = 64;
}

Kaiyu::GCOTraj::~GCOTraj() {

}

void Kaiyu::GCOTraj::Flush() {
	for (int i = 0; i < tree_nodes_.size(); i++) {
		if (tree_nodes_[i]->seg_num_ > 0) {
			tree_nodes_[i]->WritePage();
		}
	}
}

void Kaiyu::GCOTraj::InitializeDebug() {
	tree_nodes_.push_back(nullptr);
	tree_nodes_.back() = new GCOTrajNode(0);
	GCOTrajNode* root = tree_nodes_.back();
	root->depth_ = 0;
	double value[3] = { 0, 0, 0 };
	root->UpdateRange(value);
	value[0] = 10.0;
	value[1] = 10.0;
	value[2] = 10.0;
	root->UpdateRange(value);
	memcpy(boundary_min, root->min_range, sizeof(double) * 3);
	memcpy(boundary_max, root->max_range, sizeof(double) * 3);
	grid_granularity[0] = 10;
	grid_granularity[1] = 10;
	grid_granularity[2] = 10;
	for (int i = 0; i < 3; i++) {
		grid_margin[i] = (boundary_max[i] - boundary_min[i]) / grid_granularity[i];
	}
	int pos = 0;
	while (pos < tree_nodes_.size()) {
		GCOTrajNode* node = tree_nodes_[pos];
		if (node->depth_ < 3) {
			node->children_.resize(grid_granularity[node->depth_]);
			for (int i = 0; i < grid_granularity[node->depth_]; i++) {
				GCOTrajNode* child = new GCOTrajNode(tree_nodes_.size());
				node->children_[i] = child->id_;
				child->depth_ = node->depth_ + 1;
				tree_nodes_.push_back(child);
				
				for (int j = 0; j < 3; j++) {
					if (j == node->depth_) {
						child->min_range[j] = node->min_range[j] + i * (node->max_range[j] - node->min_range[j]) / grid_granularity[node->depth_];
						child->max_range[j] = child->min_range[j] + (node->max_range[j] - node->min_range[j]) / grid_granularity[node->depth_];
					}
					else {
						child->min_range[j] = node->min_range[j];
						child->max_range[j] = node->max_range[j];
					}
				}
			}
		}
		pos += 1;
	}
}
void Kaiyu::GCOTraj::Initialize(const string traj_file_name) {
	cout << "opening file " << traj_file_name << endl;
	std::ifstream ifs(traj_file_name, std::ios::in);
	unsigned int total_traj = 0;
	unsigned int point_num = 0;
	ifs >> total_traj;
	cout << total_traj << " trajectories" << endl;
	trajectories_  = vector<Trajectory*>(total_traj, nullptr);
	tree_nodes_.push_back(nullptr);
	tree_nodes_.back() = new GCOTrajNode(0);
	GCOTrajNode* root = tree_nodes_.back();
	root->depth_ = 0;
	double value[3];
	for (unsigned int i = 0; i < total_traj; i++) {
		ifs >> point_num;
		trajectories_[i] = new Trajectory(i, point_num);
		for (unsigned int point_id = 0; point_id < point_num; point_id += 1) {
//			double x, y, t;
//			ifs >> x >> y >> t;
			ifs >> value[0] >> value[1] >> value[2];
			trajectories_[i]->points_[point_id].Set(value[0], value[1], value[2]);
			root->UpdateRange(value);

		}
	}
	ifs.close();
	memcpy(boundary_min, root->min_range, sizeof(double) * 3);
	memcpy(boundary_max, root->max_range, sizeof(double) * 3);
	for (int i = 0; i < 3; i++) {
		grid_margin[i] = (boundary_max[i] - boundary_min[i]) / grid_granularity[i];
	}
	cout << "GCOTraj: " << endl;
	cout << "Granularity:";
	for (int g : grid_granularity) {
		cout << " " << g;
	}
	cout << endl;
	int pos = 0;
	while (pos < tree_nodes_.size()) {
		GCOTrajNode* node = tree_nodes_[pos];
		if (node->depth_ < 3) {
			node->children_.resize(grid_granularity[node->depth_]);
			for (int i = 0; i < grid_granularity[node->depth_]; i++) {
				GCOTrajNode* child = new GCOTrajNode(tree_nodes_.size());
				node->children_[i] = child->id_;
				child->depth_ = node->depth_ + 1;
				tree_nodes_.push_back(child);
				for (int j = 0; j < 3; j++) {
					if (j == node->depth_) {
						child->min_range[j] = node->min_range[j] + i * (node->max_range[j] - node->min_range[j]) / grid_granularity[node->depth_];
						child->max_range[j] = child->min_range[j] + (node->max_range[j] - node->min_range[j]) / grid_granularity[node->depth_];
					}
					else {
						child->min_range[j] = node->min_range[j];
						child->max_range[j] = node->max_range[j];
					}
				}
			}
		}
		pos += 1;
	}
	Segment seg;
	for (Trajectory* traj : trajectories_) {
		for (int p_id = 0; p_id < traj->points_.size() - 1; p_id++) {
			seg.Set(traj->points_[p_id], traj->points_[p_id + 1]);
			seg.original_id_ = traj->id_;
			AddSegment(seg);
		}
	}
	Flush();
}


void Kaiyu::GCOTraj::RangeQuery(double range_max[3], double range_min[3]) {
	list<int> queue;
	queue.push_back(0);
	node_access = 0;
	block_access = 0;
	leaf_access = 0;
	vector<Segment> segs;
	set<int> results;
	while (!queue.empty()) {
		int nid = queue.front();
		queue.pop_front();
		GCOTrajNode* node = tree_nodes_[nid];
		if (node->depth_ == 3) {
			leaf_access += 1;
			node_access += 1;
			for (int page_id = 0; page_id < node->pages_.size(); page_id++) {
				segs.clear();
				node->LoadPage(segs, page_id);
				block_access += 1;
				for (int seg_id = 0; seg_id < segs.size(); seg_id++) {
					if (segs[seg_id].IsOverlap(range_max, range_min)) {
						results.insert(segs[seg_id].original_id_);
					}
				}
			}
		}
		else {
			for (int child_id : node->children_) {
				GCOTrajNode* child = tree_nodes_[child_id];
				if (child->IsOverlap(range_max, range_min)) {
					queue.push_back(child_id);
				}
			}
			node_access += 1;
		}
	}
	cout << results.size() << " trajectories retrieved" << endl;
}

void Kaiyu::GCOTraj::KNNQuery(int k, double temporal_range[2], double coord[2]) {
	node_access = 0;
	leaf_access = 0;
	block_access = 0;
	Point q;
	q.Set(coord[0], coord[1], 0.0);
	GCOTrajNode* root = tree_nodes_[0];
	priority_queue<pair<double, int>, vector<pair<double, int> >, std::greater<pair<double, int> > > pqueue;
	if (root->IsTemporalOverlap(temporal_range)) {
		pqueue.emplace(root->L2Distance(q), 0);
	}
	vector<Segment> segs;
	vector<pair<double, int> > results;
	while (!pqueue.empty()) {
		if (results.size() == k && results.back().first < pqueue.top().first) {
			break;
		}
		pair<double, int> top = pqueue.top();
		pqueue.pop();
		GCOTrajNode* node = tree_nodes_[top.second];
		if (node->depth_ == 3) {
			leaf_access += 1;
			node_access += 1;
			block_access += node->pages_.size();
			for (unsigned int page_id = 0; page_id < node->pages_.size(); page_id++) {
				segs.clear();
				node->LoadPage(segs, page_id);
				for (int seg_id = 0; seg_id < segs.size(); seg_id++) {
					if (!segs[seg_id].IsTemporalOverlap(temporal_range)) {
						continue;
					}
					double distance = pDistance(q.x(), q.y(), segs[seg_id].start.x(), segs[seg_id].start.y(), segs[seg_id].end.x(), segs[seg_id].end.y());
					int original_id = segs[seg_id].original_id_;
					bool is_in = false;
					bool changed = false;
					for (int i = 0; i < results.size(); i++) {
						if (results[i].second == original_id) {
							is_in = true;
							if (results[i].first > distance) {
								results[i].first = distance;
								changed = true;
							}
							break;
						}
					}
					if (!is_in) {
						if (results.size() < k) {
							results.emplace_back(distance, original_id);
							changed = true;
						}
						else {
							if (results.back().first > distance) {
								results.back().first = distance;
								results.back().second = original_id;
								changed = true;
							}
						}
					}
					if (changed) {
						sort(results.begin(), results.end());
					}
				}
			}
		}
		else {
			node_access += 1;
			if (node->children_.empty()) {
				cout << "node has not children" << endl;
				exit(0);
			}
			for (int c_id : node->children_) {
				GCOTrajNode* child = tree_nodes_[c_id];
				if (child->IsTemporalOverlap(temporal_range)) {
					double d = child->L2Distance(q);
					pqueue.emplace(d, c_id);
				}
			}
		}
	}

}



bool Kaiyu::GCOTrajNode::IsOnBoundary(Point& p) {
	bool is_in = true;
	bool is_on_boundary = false;
	for (int i = 0; i < 3; i++) {
		if (min_range[i] <= p.value(i) && p.value(i) <= max_range[i]) {
			if (min_range[i] == p.value(i) || max_range[i] == p.value(i)) {
				is_on_boundary = true;
			}
		}
		else {
			is_in = false;
			break;
		}
	}
	return (is_in && is_on_boundary);
}







bool Kaiyu::GCOTraj::IsInCell(Segment& seg, int cell_id) {
	GCOTrajNode* node = tree_nodes_[cell_id];
	bool is_in[2] = { true, true };
	for (int dim = 0; dim < 3; dim++) {
		if (node->min_range[dim] <= seg.start.value(dim) && seg.start.value(dim) <= node->max_range[dim]) {
			//do nothing
		}
		else {
			is_in[0] = false;
		}
		if (node->min_range[dim] <= seg.end.value(dim) && seg.end.value(dim) <= node->max_range[dim]) {
			//do nothing
		}
		else {
			is_in[1] = false;
		}
	}
	if (is_in[0] && is_in[1]) {
		return true;
	}
	else {
		return false;
	}
}

void Kaiyu::GCOTraj::AddSegment(Segment& seg) {
	list<pair<int, Segment*> > seg_queue;
	double coord1[3];
	double coord2[3];
	Segment* ori_seg = new Segment(seg.start.data_addr(), seg.end.data_addr());
	ori_seg->original_id_ = seg.original_id_;
	seg_queue.emplace_back(0, ori_seg);
	while (!seg_queue.empty()) {
		pair<int, Segment*> top = seg_queue.front();
		seg_queue.pop_front();
		GCOTrajNode* node = tree_nodes_[top.first];
		if (IsInCell(*(top.second), top.first)) {
			if (node->depth_ == 3) {
				node->AddSegment(*(top.second));
			}
			else {
				int split_dim = node->depth_;
				double min_v = min(top.second->start.value(split_dim), top.second->end.value(split_dim));
				double max_v = max(top.second->start.value(split_dim), top.second->end.value(split_dim));
				for (int child_id : tree_nodes_[top.first]->children_) {
					GCOTrajNode* child = tree_nodes_[child_id];
					if (min_v <= child->min_range[split_dim] && child->max_range[split_dim] <= max_v) {
						//the segment cross the cell.
						double ratio1 = (child->min_range[split_dim] - top.second->start.value(split_dim)) / (top.second->end.value(split_dim) - top.second->start.value(split_dim));
						double ratio2 = (child->max_range[split_dim] - top.second->start.value(split_dim)) / (top.second->end.value(split_dim) - top.second->start.value(split_dim));
						for (int i = 0; i < 3; i++) {
							coord1[i] = ratio1 * (top.second->end.value(i) - top.second->start.value(i)) + top.second->start.value(i);
							coord2[i] = ratio2 * (top.second->end.value(i) - top.second->start.value(i)) + top.second->start.value(i);
						}
						Segment* subseg = new Segment(coord1, coord2);
						subseg->original_id_ = top.second->original_id_;
						seg_queue.emplace_back(child_id, subseg);
					}
					else if (min_v <= child->min_range[split_dim] && child->min_range[split_dim] < max_v && max_v < child->max_range[split_dim]) {
						double ratio = (child->min_range[split_dim] - top.second->start.value(split_dim)) / (top.second->end.value(split_dim) - top.second->start.value(split_dim));
						for (int i = 0; i < 3; i++) {
							coord1[i] = ratio * (top.second->end.value(i) - top.second->start.value(i)) + top.second->start.value(i);
						}
						Segment* subseg = nullptr;;
						if (top.second->end.value(split_dim) > top.second->start.value(split_dim)) {
							subseg = new Segment(coord1, top.second->end.data_addr());
						}
						else {
							subseg = new Segment(coord1, top.second->end.data_addr());
						}
						subseg->original_id_ = top.second->original_id_;
						seg_queue.emplace_back(child_id, subseg);
					}
					else if (child->min_range[split_dim] < min_v && min_v < child->max_range[split_dim] && child->max_range[split_dim] <= max_v) {
						double ratio = (child->max_range[split_dim] - top.second->start.value(split_dim)) / (top.second->end.value(split_dim) - top.second->start.value(split_dim));
						for (int i = 0; i < 3; i++) {
							coord1[i] = ratio * (top.second->end.value(i) - top.second->start.value(i)) + top.second->start.value(i);
						}
						Segment *subseg = nullptr;
						if (top.second->end.value(split_dim) > top.second->start.value(split_dim)) {
							subseg = new Segment(coord1, top.second->start.data_addr());
						}
						else {
							subseg = new Segment(coord1, top.second->end.data_addr());
						}
						subseg->original_id_ = top.second->original_id_;
						seg_queue.emplace_back(child_id, subseg);
					}
					else if (child->min_range[split_dim] <= min_v && max_v < child->max_range[split_dim]) {
						Segment* subseg = new Segment(top.second->start.data_addr(), top.second->end.data_addr());
						subseg->original_id_ = top.second->original_id_;
						seg_queue.emplace_back(child_id, subseg);
					}
					else if (min_v == max_v && max_v == child->max_range[split_dim] && max_v == boundary_max[split_dim]) {
						Segment* subseg = new Segment(top.second->start.data_addr(), top.second->end.data_addr());
						subseg->original_id_ = top.second->original_id_;
						seg_queue.emplace_back(child_id, subseg);
					}
					
				}
			}
		}
		delete top.second;
		//else {
		//	cout << "[Error]:" << endl;
		//	cout << "Cell boundary: ";
		//	for (int i = 0; i < 3; i++) {
		//		cout <<"(" << tree_nodes_[top.first]->min_range[i] << ", " << tree_nodes_[top.first]->max_range[i] << ") ";
		//	}
		//	cout << "Segment: [";
		//	for (int i = 0; i < 3; i++) {
		//		cout << top.second.start.value(i) << " ";
		//	}
		//	cout << "] -> [";
		//	for (int i = 0; i < 3; i++) {
		//		cout << top.second.end.value(i) << " ";
		//	}
		//	cout << endl;
		//}
	}
}


unsigned int Kaiyu::GCOTraj::LocatePointWithReference(Point& point, Point& reference) {
	int idx[3];
	unsigned int node_id;
	for (int i = 0; i < 3; i++) {
		idx[i] = (int)floor((point.value(i) - boundary_min[i]) / grid_margin[i]);
		if (point.value(i) == boundary_min[i] + idx[i] * grid_margin[i]) {
			if (reference.value(i) > point.value(i)) {
				//do nothing
			}
			else {
				idx[i] = idx[i] - 1;
			}
		}
	}
	node_id = 1 + grid_granularity[0] + grid_granularity[0] * grid_granularity[1] + idx[0] * grid_granularity[1] * grid_granularity[2] + idx[1] * grid_granularity[2] + idx[2];
	return node_id;
}

unsigned int Kaiyu::GCOTraj::LocatePoint(Point& point) {
	int idx[3];
	for (int i = 0; i < 3; i++) {
		idx[i] = (int)floor((point.value(i) - tree_nodes_[0]->min_range[i]) / grid_margin[i]);
	}
	unsigned int node_id = 1 + grid_granularity[0] + grid_granularity[0] * grid_granularity[1] + idx[0] * grid_granularity[1] * grid_granularity[2] + idx[1] * grid_granularity[2] + idx[2];
	return node_id;
}



int Kaiyu::KDTree::Query(double* query_max, double* query_min) {
	list<int> queue;
	set<int> result;
	node_access_ = 0;
	leaf_node_access_ = 0;
	if (tree_nodes_[0]->IsOverlap(query_max, query_min)) {
		dlog(cout << "root is overlapped" << endl;)
			queue.push_back(0);
	}
	while (!queue.empty()) {
		int node_id = queue.front();
		queue.pop_front();
		KDTreeNode* node = tree_nodes_[node_id];
		node_access_ += 1;
		if (node->is_leaf_) {
			leaf_node_access_ += 1;
			for (int id : node->data_) {
				Trajectory* traj = trajectories_[id];

				dlog(cout << "processing trajectory" << endl;)
					dlog(traj->Print();)

					if (traj->IsOverlap(query_max, query_min)) {
						dlog(cout << "is overlapped" << endl;)
							result.insert(traj->original_id_);
					}
					else {
						//cout << "not overlapped" << endl;
					}
			}
		}
		else {
			KDTreeNode* left = tree_nodes_[node->children[0]];
			KDTreeNode* right = tree_nodes_[node->children[1]];
			dlog(left->Print();)
				if (left->IsOverlap(query_max, query_min)) {
					//cout << "node " << node->children[0] << " is overlapped" << endl;
					queue.push_back(node->children[0]);
				}
			dlog(right->Print();)
				if (right->IsOverlap(query_max, query_min)) {
					dlog(cout << "node " << node->children[1] << " is overlapped" << endl;)
						queue.push_back(node->children[1]);
				}
		}
	}
	return result.size();
}

void Kaiyu::KDTree::PropagateSegmentNum(const unsigned int node_id) {
	KDTreeNode* node = tree_nodes_[node_id];
	node->overall_segment_ = node->segment_num_;
	int segment_num = node->segment_num_;
	while (node->parent_ != 0) {
		KDTreeNode* parent = tree_nodes_[node->parent_];
		parent->overall_segment_ += segment_num;
		node = parent;
	}
	tree_nodes_[0]->overall_segment_ += segment_num;
}

void Kaiyu::KDTree::PropagateSegmentNum() {
	for (unsigned int i = 0; i < tree_nodes_.size(); i++) {
		KDTreeNode* node = tree_nodes_[i];
		if (node->is_leaf_) {
			PropagateSegmentNum(i);
		}
	}
	propagated_ = true;
}

void Kaiyu::KDTree::SplitNode(const unsigned int node_id, const unsigned int dim, const double cut) {
	KDTreeNode* node = tree_nodes_[node_id];
	node->partition_dim_ = dim;
	node->partition_idx_ = 0;
	node->children[0] = tree_nodes_.size();
	tree_nodes_.push_back(nullptr);
	node->children[1] = tree_nodes_.size();
	tree_nodes_.push_back(nullptr);
	tree_nodes_[node->children[0]] = new KDTreeNode();
	tree_nodes_[node->children[1]] = new KDTreeNode();
	auto left = tree_nodes_[node->children[0]];
	auto right = tree_nodes_[node->children[1]];
	left->depth_ = node->depth_ + 1;
	left->parent_ = node_id;
	right->depth_ = node->depth_ + 1;
	if (right->depth_ > depth_) {
		depth_ = right->depth_;
	}
	right->parent_ = node_id;
	left->is_leaf_ = true;
	right->is_leaf_ = true;

	list<int> trajectories_2be_split;
	for (unsigned int i = 0; i < node->data_.size(); i++) {
		unsigned int traj_id = node->data_[i];
		Trajectory* traj = trajectories_[traj_id];
		int cmp_result = traj->Compare(dim, cut);
		if (cmp_result < 0) {
			left->data_.push_back(traj_id);
		}
		else if (cmp_result > 0) {
			right->data_.push_back(traj_id);
		}
		else {
			trajectories_2be_split.push_back(traj_id);
		}
	}
	list<int> sub1;
	list<int> sub2;
	for (int traj_id : trajectories_2be_split) {

		SplitTrajectory(traj_id, dim, cut, sub1, sub2);
		for (int id : sub1) {
			left->data_.push_back(id);
		}
		for (int id : sub2) {
			right->data_.push_back(id);
		}
		sub1.clear();
		sub2.clear();
	}
	node->data_.clear();
	ComputeMBR(node->children[0]);
	ComputeMBR(node->children[1]);
	ComputeMBRange(node->children[0]);
	ComputeMBRange(node->children[1]);
	ComputeSegment(node->children[0]);
	ComputeSegment(node->children[1]);
	node->segment_increase_ = left->segment_num_ + right->segment_num_ - node->segment_num_;
	right->segment_increase_ = 0;
	left->segment_increase_ = 0;
	right->overall_segment_ = 0;
	left->overall_segment_ = 0;
	node->is_leaf_ = false;
	//cout << tree_nodes_[node->children[0]]->min_range[0] << " " << tree_nodes_[node->children[0]]->min_range[1] << " " << tree_nodes_[node->children[0]]->min_range[2] << endl;
}

double Kaiyu::KDTree::ComputeMedium(const unsigned int node_id, const unsigned int dim) {
	KDTreeNode* node = tree_nodes_[node_id];
	vector<double> values;
	for (int traj_id : node->data_) {
		Trajectory* traj = trajectories_[traj_id];
		for (int i = 0; i < traj->points_.size() - 1; i++) {
			if (abs(traj->points_[i].value(dim) - traj->points_[i + 1].value(dim)) <= min_prec) {
				//if (traj->points_[i].value(dim) == traj->points_[i + 1].value(dim)) {
				continue;
			}
			values.push_back(0.5 * (traj->points_[i].value(dim) + traj->points_[i + 1].value(dim)));
		}
	}
	sort(values.begin(), values.end());
	double medium = 0;
	if (values.empty()) {
		medium = 0.5 * (node->max_range[dim] + node->min_range[dim]);
	}
	else {
		medium = values[values.size() / 2];
	}
	//double medium = values[values.size() / 2];
	return medium;
}

void Kaiyu::KDTree::GetChildSegNumAhead(const unsigned int node_id, const unsigned int dim, const unsigned int cut_idx, int* child_seg_num) {
	KDTreeNode* node = tree_nodes_[node_id];
	double cut = 0;
	if (node->idx_range_max[dim] - node->idx_range_min[dim] == 1) {
		cut = ComputeMedium(node_id, dim);
	}
	else {
		cut = cuts_[dim][cut_idx];
	}
	child_seg_num[0] = 0;
	child_seg_num[1] = 0;
	for (int traj_id : node->data_) {
		Trajectory* traj = trajectories_[traj_id];
		int traj_side = 0;
		if (traj->points_[0].value(dim) < cut) {
			traj_side = 1;
		}
		else if (traj->points_[0].value(dim) > cut) {
			traj_side = 2;
		}
		int seg_num = 0;
		for (unsigned int pid = 1; pid < traj->points_.size(); pid++) {
			seg_num += 1;
			if (traj_side == 0) {
				if (traj->points_[pid].value(dim) < cut) {
					traj_side = 1;
				}
				if (traj->points_[pid].value(dim) > cut) {
					traj_side = 2;
				}
			}
			else {
				int point_side = 0;
				if (traj->points_[pid].value(dim) < cut) {
					point_side = 1;
				}
				if (traj->points_[pid].value(dim) > cut) {
					point_side = 2;
				}
				if (traj_side + point_side == 3) {
					if (traj->points_[pid - 1].value(dim) == cut) {
						child_seg_num[traj_side - 1] += seg_num - 1;
						traj_side = point_side;
						seg_num = 1;

					}
					else {
						child_seg_num[traj_side - 1] += seg_num;
						seg_num = 1;
						traj_side = point_side;
					}
				}
			}
		}
	}
}

int Kaiyu::KDTree::GetSegmentIncreaseAhead(const unsigned int node_id, const unsigned int dim, const unsigned int cut_idx) {
	int increase = 0;
	KDTreeNode* node = tree_nodes_[node_id];
	double cut = 0;
	if (node->idx_range_max[dim] - node->idx_range_min[dim] == 1) {
		cut = ComputeMedium(node_id, dim);
	}
	else {
		cut = cuts_[dim][cut_idx];
	}
	for (int traj_id : node->data_) {
		Trajectory* traj = trajectories_[traj_id];
		for (unsigned pid = 0; pid < traj->points_.size() - 1; pid++) {
			double min_value = min(traj->points_[pid].value(dim), traj->points_[pid + 1].value(dim));
			double max_value = max(traj->points_[pid].value(dim), traj->points_[pid + 1].value(dim));
			if (min_value < cut && cut < max_value) {
				increase += 1;
			}
		}
	}
	return increase;
}

int Kaiyu::KDTree::SplitNode(const unsigned int node_id, const unsigned int dim, const unsigned int cut_idx) {
	//cout << "cut idx: " << cut_idx << endl;
	bool debug = false;
	KDTreeNode* node = tree_nodes_[node_id];
	if (node_id == 826699) {
		cout << "into debugging" << endl;
		debug = true;
	}
	int split_by_predefined_cuts = 0;
	int split_dim = dim;
	double cut = 0;
	if (node->idx_range_max[dim] - node->idx_range_min[dim] == 1) {
		if (!IsSplitable(node_id, dim)) {
			if (debug) {
				cout << "dimension " << dim << " can not be splitted" << endl;
			}
			for (int d = 0; d < 3; d++) {
				if (d == dim)continue;
				if (IsSplitable(node_id, d)) {
					split_dim = d;
					if (debug) {
						cout << "Split in the " << d << " dimension instead" << endl;
					}
					//
					//getchar();
					break;
				}
			}
		}
		cut = ComputeMedium(node_id, split_dim);
		if (debug) {
			cout << "[debug] coords:" << cut - node->min_range[split_dim] << " " << node->max_range[split_dim] - node->min_range[split_dim] << endl;
			PrintTrajectoryDiff(node_id, split_dim, cut);
		}
		SplitNode(node_id, split_dim, cut);


	}
	else {
		cut = cuts_[dim][cut_idx];
		SplitNode(node_id, dim, cut);
		split_by_predefined_cuts = 1;
		node->partition_idx_ = cut_idx;
	}
	auto left = tree_nodes_[node->children[0]];
	auto right = tree_nodes_[node->children[1]];
	left->depth_ = node->depth_ + 1;
	left->parent_ = node_id;
	right->depth_ = node->depth_ + 1;
	if (right->depth_ > depth_) {
		depth_ = right->depth_;
	}
	right->parent_ = node_id;
	left->is_leaf_ = true;
	right->is_leaf_ = true;
	if (split_by_predefined_cuts == 1) {
		right->is_splitted_by_idx_ = true;
		left->is_splitted_by_idx_ = true;
		/*for (unsigned int i = 0; i < 3; i++) {
			if (i == dim) {
				left->idx_range_min[i] = node->idx_range_min[i];
				left->idx_range_max[i] = cut_idx;
				right->idx_range_min[i] = cut_idx;
				right->idx_range_max[i] = node->idx_range_max[i];
			}
			else {
				left->idx_range_min[i] = node->idx_range_min[i];
				left->idx_range_max[i] = node->idx_range_max[i];
				right->idx_range_min[i] = node->idx_range_min[i];
				right->idx_range_max[i] = node->idx_range_max[i];
			}
		}*/
	}
	else {
		right->is_splitted_by_idx_ = false;
		left->is_splitted_by_idx_ = false;
		/*for (unsigned int i = 0; i < 3; i++) {
			left->idx_range_min[i] = node->idx_range_min[i];
			left->idx_range_max[i] = node->idx_range_max[i];
			right->idx_range_min[i] = node->idx_range_min[i];
			right->idx_range_max[i] = node->idx_range_max[i];
		}*/
	}
	if (left->data_.size() == 0 || right->data_.size() == 0) {
		cout << "One of the children has zero segments" << endl;
		cout << "left: " << left->max_range[0] << " " << left->max_range[1] << " " << left->min_range[0] << " " << left->min_range[1] << endl;
		cout << "right: " << right->max_range[0] << " " << right->max_range[1] << " " << right->min_range[0] << " " << right->min_range[1] << endl;
		cout << "cut: " << cut << endl;
		cout << "segment num: " << left->segment_num_ << " " << right->segment_num_ << endl;

		exit(0);
	}
	return split_by_predefined_cuts;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

void Kaiyu::KDTree::PrintTrajectoryDiff(const unsigned int node_id, const unsigned int dim, const double cut) {
	Kaiyu::KDTreeNode* node = tree_nodes_[node_id];
	cout << node->data_.size() << " trajectories: " << endl;
	for (int i = 0; i < node->data_.size(); i++) {
		int traj_id = node->data_[i];
		Kaiyu::Trajectory* traj = trajectories_[traj_id];
		cout << "Trajectory " << i << ":";
		for (Kaiyu::Point& p : traj->points_) {
			cout << p.value(dim) - cut << " ";
		}
		cout << endl;
	}
}

void Kaiyu::KDTree::PrintTrajectory(const unsigned int node_id) {
	Kaiyu::KDTreeNode* node = tree_nodes_[node_id];
	cout << node->data_.size() << " trajectories: " << endl;
	for (int i = 0; i < node->data_.size(); i++) {
		int traj_id = node->data_[i];
		Kaiyu::Trajectory* traj = trajectories_[traj_id];
		cout << "Trajectory " << i << ":";
		for (Kaiyu::Point& p : traj->points_) {
			cout << " (" << p.value(0) << ", " << p.value(1) << ", " << p.value(2) << ")";
		}
		cout << endl;
	}
}

Kaiyu::TBTreeNode::TBTreeNode() {
	min_range[0] = DBL_MAX;
	min_range[1] = DBL_MAX;
	min_range[2] = DBL_MAX;
	max_range[0] = -DBL_MAX;
	max_range[1] = -DBL_MAX;
	max_range[2] = -DBL_MAX;
}

void Kaiyu::TBTreeNode::UpdateRangeFromNode(TBTreeNode* node) {
	for (int dim = 0; dim < 3; dim++) {
		if (node->max_range[dim] > max_range[dim]) {
			max_range[dim] = node->max_range[dim];
		}
		if (node->min_range[dim] < min_range[dim]) {
			min_range[dim] = node->min_range[dim];
		}
	}
}

void Kaiyu::TBTreeNode::UpdateRange(Point& p) {
	if (p.x() > max_range[0]) {
		max_range[0] = p.x();
	}
	if (p.x() < min_range[0]) {
		min_range[0] = p.x();
	}
	if (p.y() > max_range[1]) {
		max_range[1] = p.y();
	}
	if (p.y() < min_range[1]) {
		min_range[1] = p.y();
	}
	if (p.t() > max_range[2]) {
		max_range[2] = p.t();
	}
	if (p.t() < min_range[2]) {
		min_range[2] = p.t();
	}
}

void TBTreeRangeQuery(Kaiyu::TBTree* tree, double query_max[3], double query_min[3]) {
	tree->RangeQuery(query_max, query_min);
}

void GCORangeQuery(Kaiyu::GCOTraj* gco, double query_max[3], double query_min[3]) {
	gco->RangeQuery(query_max, query_min);
}


Kaiyu::TBTree::~TBTree() {
	for (int i = 0; i < original_trajectories_.size(); i++) {
		delete original_trajectories_[i];
	}
	for (int i = 0; i < trajectories_.size(); i++) {
		delete trajectories_[i];
	}
	for (int i = 0; i < tree_nodes_.size(); i++) {
		delete tree_nodes_[i];
	}
}




void Kaiyu::TrajStore::KNNQuery(int k, double temporal_interval[2], double coord[2]) {
	node_access = 0;
	leaf_access = 0;
	block_acces = 0;
	Point q;
	q.Set(coord[0], coord[1], 0.0);
	TrajStoreNode* root = tree_nodes_[0];
	priority_queue<pair<double, int>, vector<pair<double, int> >, std::greater<pair<double, int> > > pqueue;
	if (root->IsTemporalOverlap(temporal_interval)) {
		pqueue.emplace(root->L2Distance(q), 0);
	}
	vector<pair<double, int> > results;
	while (!pqueue.empty()) {
		if (results.size() == k && results.back().first < pqueue.top().first) {
			break;
		}
		pair<double, int> top = pqueue.top();
		pqueue.pop();
		TrajStoreNode* node = tree_nodes_[top.second];
		if (node->is_leaf_) {
			leaf_access += 1;
			node_access += 1;
			block_acces += (int)ceil(1.0 * node->seg_num_ * (sizeof(unsigned int) + sizeof(double) * 6) / page_size_);
			Segment seg;
			unsigned int pos = 0;
			for (int i = 0; i < node->seg_num_; i++) {
				unsigned int l = seg.Deserialization(node->content + pos);
				pos += l;
				if (!seg.IsTemporalOverlap(temporal_interval)) {
					continue;
				}
				double distance = pDistance(q.x(), q.y(), seg.start.x(), seg.start.y(), seg.end.x(), seg.end.y());
				int original_id = seg.original_id_;
				bool is_in = false;
				bool changed = false;
				for (int i = 0; i < results.size(); i++) {
					if (results[i].second == original_id) {
						is_in = true;
						if (results[i].first > distance) {
							results[i].first = distance;
							changed = true;
						}
						break;
					}
				}
				if (!is_in) {
					if (results.size() < k) {
						results.emplace_back(distance, original_id);
						changed = true;
					}
					else {
						if (results.back().first > distance) {
							results.back().first = distance;
							results.back().second = original_id;
							changed = true;
						}
					}
				}
				if (changed) {
					sort(results.begin(), results.end());
				}
			}
		}
		else {
			node_access += 1;
			if (node->children_ == nullptr) {
				cout << "node has not children" << endl;
				exit(0);
			}
			for (int c_id : *(node->children_)) {
				TrajStoreNode* child = tree_nodes_[c_id];
				if (child->IsTemporalOverlap(temporal_interval)) {
					double d = child->L2Distance(q);
					pqueue.emplace(d, c_id);
				}
			}
		}
	}
}

void Kaiyu::TrajStore::RangeQuery(double query_max[3], double query_min[3]) {
	node_access = 0;
	leaf_access = 0;
	block_acces = 0;
	list<int> node_queue;
	TrajStoreNode* root = tree_nodes_[0];
	if (root->IsOverlap(query_max, query_min)) {
		node_queue.push_back(0);
	}
	set<int> results;
	while (!node_queue.empty()) {
		int node_id = node_queue.front();
		node_queue.pop_front();
		TrajStoreNode* node = tree_nodes_[node_id];
		if (node->is_leaf_) {
			Segment seg;
			unsigned int pos = 0;
			for (int seg_id = 0; seg_id < node->seg_num_; seg_id++) {
				unsigned int l = seg.Deserialization(node->content + pos);
				pos += l;
				if (seg.IsOverlap(query_max, query_min)) {
					results.insert(seg.original_id_);
				}
			}
			node_access += 1;
			leaf_access += 1;
			block_acces += (int)ceil(1.0 * node->seg_num_ * (sizeof(unsigned int) + sizeof(double) * 6) / page_size_);
		}
		else {
			if (node->children_ == nullptr) {
				cout << "node has no children" << endl;
				exit(0);
			}
			for (int c_id : *(node->children_)) {
				TrajStoreNode* child = tree_nodes_[c_id];
				if (child->IsOverlap(query_max, query_min)) {
					node_queue.push_back(c_id);
				}
			}
			node_access += 1;
		}
	}
	cout << results.size() << " trajectories retrieved" << endl;
	/*for (int i : results) {
		cout << i << " ";
	}
	cout << endl;*/
}

void Kaiyu::TBTree::RangeQuery(double range_max[3], double range_min[3]) {
	node_access = 0;
	leaf_access = 0;
	list<int> node_queue;
	TBTreeNode* root = tree_nodes_[root_];
	if (root->IsOverlap(range_max, range_min)) {
		node_queue.push_back(root_);
	}
	set<int> results;
	while (!node_queue.empty()) {
		int node_id = node_queue.front();
		node_queue.pop_front();
		TBTreeNode* node = tree_nodes_[node_id];
		if (node->is_leaf_) {
			Trajectory* traj = trajectories_[node->traj_pos_];
			cout << "(" << traj->points_[0].x() << " " << traj->points_[0].y() << ") (" << traj->points_[1].x() << " " << traj->points_[1].y() << ")" << endl;
			if (traj->IsOverlap(range_max, range_min)) {
				cout << "overlap" << endl;
				results.insert(traj->original_id_);
			}
			leaf_access += 1;
			node_access += 1;
		}
		else {
			for (int child_id : node->children_) {
				TBTreeNode* child = tree_nodes_[child_id];
				if (child->IsOverlap(range_max, range_min)) {
					node_queue.push_back(child_id);
				}
			}
			node_access += 1;
		}
	}
	cout << results.size() << " trajectories retrieved" << endl;
	for (int id : results) {
		cout << id << " ";
	}
	cout << endl;
}

void TSRangeQuery(Kaiyu::TrajStore* ts, double query_max[3], double query_min[3]) {
	ts->RangeQuery(query_max, query_min);
}

void TSKNNQuery(Kaiyu::TrajStore* ts, int k, double temporal_interval[2], double coord[2]) {
	ts->KNNQuery(k, temporal_interval, coord);
}

void GCOKNNQuery(Kaiyu::GCOTraj* gco, int k, double temporal_range[2], double coord[2]) {
	gco->KNNQuery(k, temporal_range, coord);
}

int GetTSNodeAccess(Kaiyu::TrajStore* ts) {
	return ts->node_access;
}

int GetTSLeafAccess(Kaiyu::TrajStore* ts) {
	return ts->leaf_access;
}

int GetTSBlockAccess(Kaiyu::TrajStore* ts) {
	return ts->block_acces;
}

int GetGCONodeAccess(Kaiyu::GCOTraj* gco) {
	return gco->node_access;
}

int GetGCOLeafAccess(Kaiyu::GCOTraj* gco) {
	return gco->leaf_access;
}

void SetGCOGranularity(Kaiyu::GCOTraj* gco, int lv1, int lv2, int lv3) {
	gco->grid_granularity[0] = lv1;
	gco->grid_granularity[1] = lv2;
	gco->grid_granularity[2] = lv3;
}

void SetGCOSegNumPerPage(Kaiyu::GCOTraj* gco, int num) {
	GCOTrajNode_Segnum = num;
}

int GetGCOBlockAccess(Kaiyu::GCOTraj* gco) {
	return gco->block_access;
}

void SetTSStopIncreaseRatio(Kaiyu::TrajStore* ts, double ratio) {
	ts->stop_increase_ratio = ratio;
}

void SetTSPageSize(Kaiyu::TrajStore* ts, unsigned int page_size) {
	ts->B = page_size;
}

void SetTSQueryPara(Kaiyu::TrajStore* ts, double qx, double qy) {
	ts->query_spatial_domain[0] = qx;
	ts->query_spatial_domain[1] = qy;
}
void SetTSEpsilon(Kaiyu::TrajStore* ts, double epsilon) {
	ts->epsilon = epsilon;
}
void SetTSFanout(Kaiyu::TrajStore* ts, int fanout) {
	ts->B = fanout;
}

int GetTBTreeNodeAccess(Kaiyu::TBTree* tree) {
	return tree->node_access;
}
int GetTBTreeLeafAccess(Kaiyu::TBTree* tree) {
	return tree->leaf_access;
}


double pDistance(double x, double y, double x1, double y1, double x2, double y2) {
	double a = x - x1;
	double b = y - y1;
	double c = x2 - x1;
	double d = y2 - y1;

	double dot = a * c + b * d;
	double len_sq = c * c + d * d;
	double param = -1;
	if (len_sq != 0) {
		param = dot / len_sq;
	}
	double xx, yy;
	if (param < 0) {
		xx = x1;
		yy = y1;
	}
	else if (param > 1) {
		xx = x2;
		yy = y2;
	}
	else {
		xx = x1 + param * c;
		yy = y1 + param * d;
	}
	double dx = x - xx;
	double dy = y - yy;
	return sqrt(dx * dx + dy * dy);
}

double* Kaiyu::Point::data_addr() {
	return values_;
}

double Kaiyu::Trajectory::DistanceToPoint(Point& p) {
	double min_dist = DBL_MAX;
	for (int i = 0; i < points_.size() - 1; i++) {
		double d = pDistance(p.x(), p.y(), points_[i].x(), points_[i].y(), points_[i + 1].x(), points_[i + 1].y());
		if (d < min_dist) {
			min_dist = d;
		}
	}
	return min_dist;
}



unsigned int Kaiyu::Trajectory::Serielization(char* content) {
	unsigned int pos = 0;
	memcpy(content + pos, &original_id_, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(content + pos, &id_, sizeof(id_));
	pos += sizeof(id_);
	unsigned int size = points_.size();
	memcpy(content + pos, &size, sizeof(size));
	pos += sizeof(size);
	for (int i = 0; i < points_.size(); i++) {
		memcpy(content + pos, points_[i].data_addr(), sizeof(double) * 3);
		pos += sizeof(double) * 3;
	}
	return pos;
}

unsigned int Kaiyu::Trajectory::ByteSize()const {
	unsigned int size = (points_.size() - 1) * (sizeof(unsigned int) + sizeof(double) * 6);
	return size;
}

unsigned int Kaiyu::Trajectory::Deserielization(char* content) {
	unsigned int pos = 0;
	memcpy(&original_id_, content + pos, sizeof(original_id_));
	pos += sizeof(original_id_);
	memcpy(&id_, content + pos, sizeof(id_));
	pos += sizeof(id_);
	unsigned int size;
	memcpy(&size, content + pos, sizeof(size));
	pos += sizeof(size);
	points_ = vector<Point>(size);
	for (int i = 0; i < size; i++) {
		memcpy(points_[i].data_addr(), content + pos, sizeof(double) * 3);
		pos += sizeof(double) * 3;
	}
	return pos;
}

double Kaiyu::KDTreeNode::L2Distance(Point& p) {
	if (min_range[0] < p.x() && p.x() < max_range[0]) {
		if (p.y() > max_range[1]) {
			return p.y() - max_range[1];
		}
		if (p.y() < min_range[1]) {
			return min_range[1] - p.y();
		}
		return 0;
	}
	else {
		if (min_range[1] < p.y() && p.y() < max_range[1]) {
			if (p.x() > max_range[0]) {
				return p.x() - max_range[0];
			}
			else {
				return min_range[0] - p.x();
			}
		}
		double dist = 0;
		double anchor[2];
		if (p.x() > max_range[0]) {
			anchor[0] = max_range[0];
		}
		else {
			anchor[0] = min_range[0];
		}
		if (p.y() > max_range[1]) {
			anchor[1] = max_range[1];
		}
		else {
			anchor[1] = min_range[1];
		}
		dist = sqrt((anchor[0] - p.x()) * (anchor[0] - p.x()) + (anchor[1] - p.y()) * (anchor[1] - p.y()));
		return dist;
	}
}
double Kaiyu::TBTreeNode::L2Distance(Point& p) {
	if (min_range[0] < p.x() && p.x() < max_range[0]) {
		if (p.y() > max_range[1]) {
			return p.y() - max_range[1];
		}
		if (p.y() < min_range[1]) {
			return min_range[1] - p.y();
		}
		return 0;
	}
	else {
		if (min_range[1] < p.y() && p.y() < max_range[1]) {
			if (p.x() > max_range[0]) {
				return p.x() - max_range[0];
			}
			else {
				return min_range[0] - p.x();
			}
		}
		double dist = 0;
		double anchor[2];
		if (p.x() > max_range[0]) {
			anchor[0] = max_range[0];
		}
		else {
			anchor[0] = min_range[0];
		}
		if (p.y() > max_range[1]) {
			anchor[1] = max_range[1];
		}
		else {
			anchor[1] = min_range[1];
		}
		dist = sqrt((anchor[0] - p.x()) * (anchor[0] - p.x()) + (anchor[1] - p.y()) * (anchor[1] - p.y()));
		return dist;
	}
}

void TBTreeKNNQuery(Kaiyu::TBTree* tree, int k, double temporal_range[2], double coord[2]) {
	tree->KNNQuery(k, temporal_range, coord);
}

void Kaiyu::TBTree::KNNQuery(int k, double temporal_range[2], double coord[2]) {
	node_access = 0;
	leaf_access = 0;
	TBTreeNode* root = tree_nodes_[root_];
	Point q;
	q.Set(coord[0], coord[1], 0.0);
	priority_queue<pair<double, int>, vector<pair<double, int> >, std::greater<pair<double, int> > > pqueue;
	if (root->IsTemporalOverlap(temporal_range)) {
		double dist = root->L2Distance(q);
		pqueue.emplace(dist, root_);
	}
	vector<pair<double, int> > results;
	while (!pqueue.empty()) {
		if (results.size() == k && results.back().first < pqueue.top().first) {
			break;
		}
		pair<double, int> top = pqueue.top();
		pqueue.pop();
		TBTreeNode* node = tree_nodes_[top.second];
		if (node->is_leaf_) {
			leaf_access += 1;
			node_access += 1;
			Trajectory* traj = trajectories_[node->traj_pos_];
			if (!traj->IsTemporalOverlap(temporal_range)) {
				continue;
			}
			double distance = traj->DistanceToPoint(q);
			int original_id = traj->original_id_;
			bool is_in = false;
			bool changed = false;
			for (int i = 0; i < results.size(); i++) {
				if (results[i].second == original_id) {
					is_in = true;
					if (results[i].first > distance) {
						results[i].first = distance;
						changed = true;
					}
					break;
				}
			}
			if (!is_in) {
				if (results.size() < k) {
					results.emplace_back(distance, original_id);
					changed = true;
				}
				else {
					if (results.back().first > distance) {
						results.back().first = distance;
						results.back().second = original_id;
						changed = true;
					}
				}
			}
			if (changed) {
				sort(results.begin(), results.end());
			}
		}
		else {
			node_access += 1;
			for (int child_id : node->children_) {
				TBTreeNode* child = tree_nodes_[child_id];
				if (child->IsTemporalOverlap(temporal_range)) {
					double d = child->L2Distance(q);
					pqueue.emplace(d, child_id);
				}
			}
		}
	}
}

void Kaiyu::TBTree::Initialize(const string traj_file_name) {
	cout << "opening file " << traj_file_name << endl;
	std::ifstream ifs(traj_file_name, std::ios::in);
	unsigned int total_traj = 0;
	unsigned int point_num = 0;
	ifs >> total_traj;
	cout << total_traj << " trajectories" << endl;
	original_trajectories_ = vector<Trajectory*>(total_traj, nullptr);
	for (unsigned int i = 0; i < total_traj; i++) {
		ifs >> point_num;
		original_trajectories_[i] = new Trajectory(i, point_num);
		for (unsigned int point_id = 0; point_id < point_num; point_id += 1) {
			double x, y, t;
			ifs >> x >> y >> t;
			original_trajectories_[i]->points_[point_id].Set(x, y, t);
		}
	}
	ifs.close();
	cout << "trajectories loaded" << endl;
	//cout << "traj_num_per_leaf: " << traj_num_per_leaf << endl;
	//cout << "child_num_per_node: " << child_num_per_node << endl;
	int sub_traj_id = 0;
	for (int traj_id = 0; traj_id < total_traj; traj_id++) {
		Trajectory* traj = original_trajectories_[traj_id];
		int length = traj->points_.size();
		int start_pos = 0;
		while (start_pos < traj->points_.size() - 1) {
			//cout << "beggining length: " << length << endl;
			int sub_traj_length = (int)min(length, traj_num_per_leaf);
			Trajectory* sub_traj = new Trajectory(traj_id, sub_traj_length);
			TBTreeNode* node = new TBTreeNode();
			sub_traj->id_ = sub_traj_id;
			for (int point_id = 0; point_id < sub_traj_length; point_id++) {
				sub_traj->points_[point_id].Set(traj->points_[start_pos + point_id]);
				node->UpdateRange(sub_traj->points_[point_id]);
			}
			node->traj_pos_ = sub_traj_id;
			node->id_ = sub_traj_id;
			node->is_leaf_ = true;
			length = length - sub_traj_length + 1;
			start_pos = start_pos + sub_traj_length - 1;
			trajectories_.push_back(sub_traj);
			tree_nodes_.push_back(node);
			sub_traj_id += 1;
			//cout << "ending length: " << length << endl;
		}
	}
	vector<TBTreeNode*> nodes_for_packing;
	nodes_for_packing.assign(tree_nodes_.begin(), tree_nodes_.end());
	vector<TBTreeNode*> nodes_next_level;
	int level = 1;
	int node_id = sub_traj_id;
	while (!nodes_for_packing.empty()) {
		//cout << "nodes_for_packing.size(): " << nodes_for_packing.size() << endl;
		int cmp_dim = level % 2;
		sort(nodes_for_packing.begin(), nodes_for_packing.end(), [cmp_dim](const TBTreeNode* node1, const TBTreeNode* node2) {return 0.5 * (node1->max_range[cmp_dim] + node1->min_range[cmp_dim]) < 0.5 * (node2->max_range[cmp_dim] + node2->min_range[cmp_dim]); });
		int start_pos = 0;
		while (start_pos < nodes_for_packing.size()) {
			TBTreeNode* node = new TBTreeNode();
			node->traj_pos_ = -1;
			node->id_ = node_id;
			node->is_leaf_ = false;
			int node_entry_num = child_num_per_node;
			if (start_pos + child_num_per_node > nodes_for_packing.size()) {
				node_entry_num = nodes_for_packing.size() - start_pos;
			}
			node->children_ = vector<int>(node_entry_num);
			for (int i = 0; i < node_entry_num; i++) {
				node->UpdateRangeFromNode(nodes_for_packing[start_pos + i]);
				node->children_[i] = nodes_for_packing[start_pos + i]->id_;
			}
			tree_nodes_.push_back(node);
			nodes_next_level.push_back(node);
			start_pos += node_entry_num;
			node_id += 1;
		}
		level += 1;
		nodes_for_packing.clear();
		if (nodes_next_level.size() > 1) {
			nodes_for_packing.assign(nodes_next_level.begin(), nodes_next_level.end());
			nodes_next_level.clear();
		}
	}
	/*for (int i = 0; i < tree_nodes_.size(); i++) {
		cout << "node " << i <<" "<<tree_nodes_[i]->id_<< endl;
		for (int dim = 0; dim < 3; dim++) {
			cout << "(" << tree_nodes_[i]->min_range[dim] << ", " << tree_nodes_[i]->max_range[dim] << ") ";
		}
		cout << endl;
		if (tree_nodes_[i]->is_leaf_) {
			cout << "traj_pos: " << tree_nodes_[i]->traj_pos_ << endl;
		}
		else {
			cout << "children:";
			for (int c : tree_nodes_[i]->children_) {
				cout << " " << c;
			}
			cout << endl;
		}
	}
	for (int i = 0; i < trajectories_.size(); i++) {
		cout << "trajctory " << i << " " << trajectories_[i]->id_ << endl;
		for (Point& p : trajectories_[i]->points_) {
			cout << "(";
			for (int dim = 0; dim < 3; dim++) {
				cout << p.value(dim) << " ";
			}
			cout << ") ";
		}
		cout << endl;
	}*/
	root_ = tree_nodes_.size() - 1;
}

Kaiyu::KDTree* ConstructTree() {
	Kaiyu::KDTree* tree = new Kaiyu::KDTree();
	return tree;
}

Kaiyu::GCOTraj* ConstructGCOTraj() {
	Kaiyu::GCOTraj* gco = new Kaiyu::GCOTraj();
	return gco;
}



Kaiyu::TBTree* ConstructTBTree(int traj_num, int child_num) {
	Kaiyu::TBTree* tree = new Kaiyu::TBTree();
	tree->traj_num_per_leaf = traj_num;
	tree->child_num_per_node = child_num;
	return tree;
}

Kaiyu::TrajStore* ConstructTrajStore() {
	Kaiyu::TrajStore* ts = new Kaiyu::TrajStore();
	ts->query_spatial_domain[0] = 100.0;
	ts->query_spatial_domain[1] = 100.0;
	ts->epsilon = 0;
	ts->B = 50;
	return ts;
}

void InitializeTBTree(Kaiyu::TBTree* tree, char* traj_file_name) {
	string filename(traj_file_name);
	tree->Initialize(filename);
	cout << tree->original_trajectories_.size() << " trajectories are loaded." << endl;
}

void InitializeTrajStore(Kaiyu::TrajStore* ts, char* traj_file_name) {
	string filename(traj_file_name);
	ts->Initialize(filename);
}

void InitializeGCOTraj(Kaiyu::GCOTraj* gco, char* traj_file_name) {
	string filename(traj_file_name);
	gco->Initialize(filename);
}



void InitializeTree(Kaiyu::KDTree* tree) {
	ifstream ifs("config.txt", std::ios::in);
	string traj_file_str;
	string cut_file_str;
	//std::getline(ifs, traj_file_str);
	//std::getline(ifs, cut_file_str);
	ifs >> traj_file_str;
	ifs >> cut_file_str;
	ifs.close();
	cout << traj_file_str << endl;
	cout << cut_file_str << endl;
	tree->Initialize(traj_file_str, cut_file_str);
	cout << tree->trajectories_.size() << endl;
}

int SplitNode(Kaiyu::KDTree* tree, unsigned int node_id, unsigned int dim, unsigned int cut_idx) {
	return tree->SplitNode(node_id, dim, cut_idx);
}

int GetDepth(Kaiyu::KDTree* tree, unsigned int node_id) {
	return tree->tree_nodes_[node_id]->depth_;
}

int GetSegmentNum(Kaiyu::KDTree* tree, unsigned int node_id) {
	return tree->tree_nodes_[node_id]->segment_num_;
}

void GetChildren(Kaiyu::KDTree* tree, unsigned int node_id, int* status) {
	if (tree->tree_nodes_[node_id]->is_leaf_) {
		status[0] = 0;
	}
	else {
		status[0] = 1;
		status[1] = tree->tree_nodes_[node_id]->children[0];
		status[2] = tree->tree_nodes_[node_id]->children[1];
	}
}

int GetImmediateSegmentIncrease(Kaiyu::KDTree* tree, unsigned int node_id) {
	return tree->tree_nodes_[node_id]->segment_increase_;
}

void GetBoundingBox(Kaiyu::KDTree* tree, unsigned int node_id, double* bounding_max, double* bounding_min) {
	Kaiyu::KDTreeNode* node = tree->tree_nodes_[node_id];
	bounding_min[0] = node->min_range[0];
	bounding_min[1] = node->min_range[1];
	bounding_min[2] = node->min_range[2];
	bounding_max[0] = node->max_range[0];
	bounding_max[1] = node->max_range[1];
	bounding_max[2] = node->max_range[2];
}

void GetSubrange(Kaiyu::KDTree* tree, unsigned int node_id, int* range_max, int* range_min) {
	Kaiyu::KDTreeNode* node = tree->tree_nodes_[node_id];
	range_min[0] = node->idx_range_min[0];
	range_min[1] = node->idx_range_min[1];
	range_min[2] = node->idx_range_min[2];
	range_max[0] = node->idx_range_max[0];
	range_max[1] = node->idx_range_max[1];
	range_max[2] = node->idx_range_max[2];
}

int GetNodeAccess(Kaiyu::KDTree* tree) {
	return tree->node_access_;
}

int GetLeafNodeAccess(Kaiyu::KDTree* tree) {
	return tree->leaf_node_access_;
}

int GetOverallSegmentIncrease(Kaiyu::KDTree* tree, unsigned int node_id) {
	if (!tree->propagated_) {
		cout << "The segment number has not been propagated yet." << endl;
		return 0;
	}
	return tree->tree_nodes_[node_id]->overall_segment_ - tree->tree_nodes_[node_id]->segment_num_;
}

void Propagate(Kaiyu::KDTree* tree) {
	tree->PropagateSegmentNum();
}

int GetPartitionDim(Kaiyu::KDTree* tree, unsigned int node_id) {
	return tree->tree_nodes_[node_id]->partition_dim_;
}

int GetPartitionIdx(Kaiyu::KDTree* tree, unsigned int node_id) {
	return tree->tree_nodes_[node_id]->partition_idx_;
}

int IsLeaf(Kaiyu::KDTree* tree, unsigned int node_id) {
	if (tree->tree_nodes_[node_id]->is_leaf_) {
		return 1;
	}
	else {
		return 0;
	}
}

int Query(Kaiyu::KDTree* tree, double* query_max, double* query_min) {
	int result = tree->Query(query_max, query_min);
	return result;
}

void KNNQuery(Kaiyu::KDTree* tree, int k, double temporal_range[2], double coord[2]) {
	tree->KNNQuery(k, temporal_range, coord);
}

void GetChildSegmentNumAhead(Kaiyu::KDTree* tree, unsigned int node_id, unsigned int dim, unsigned int cut_idx, int* child_seg_num) {
	tree->GetChildSegNumAhead(node_id, dim, cut_idx, child_seg_num);
}
int GetImmediateSegmentIncreaseAhead(Kaiyu::KDTree* tree, unsigned int node_id, unsigned int dim, unsigned int cut_idx) {
	int increase = tree->GetSegmentIncreaseAhead(node_id, dim, cut_idx);
	return increase;
}

void PrintTrajectory(Kaiyu::KDTree* tree, const unsigned int node_id) {
	tree->PrintTrajectory(node_id);
}




