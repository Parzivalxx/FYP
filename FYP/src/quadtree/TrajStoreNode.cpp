#include "TrajStoreNode.h"

FYP::TrajStoreNode::TrajStoreNode() {
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

bool FYP::TrajStoreNode::isTemporalOverlap(const double temporal_interval[2]) {
	if (max_range[2] < temporal_interval[0] || min_range[2] > temporal_interval[1]) {
		return false;
	}
	return true;
}

bool FYP::TrajStoreNode::isOverlap(const double query_max[3], const double query_min[3]) {
	for (int i = 0; i < 3; i++) {
		if (max_range[i] < query_min[i] || min_range[i] > query_max[i]) {
			return false;
		}
	}
	return true;
}

void FYP::TrajStoreNode::updateRange(double x, double y, double t) {
	max_range[0] = std::max(max_range[0], x);
	min_range[0] = std::min(min_range[0], x);
	max_range[1] = std::max(max_range[1], y);
	min_range[1] = std::min(min_range[1], y);
	max_range[2] = std::max(max_range[2], t);
	min_range[2] = std::min(min_range[2], t);
}

double FYP::TrajStoreNode::getL2Distance(Point& p) {
	if (min_range[0] < p.getX() && p.getX() < max_range[0]) {
		if (p.getY() > max_range[1]) {
			return p.getY() - max_range[1];
		}
		if (p.getY() < min_range[1]) {
			return min_range[1] - p.getY();
		}
		return 0;
	}
	else {
		if (min_range[1] < p.getY() && p.getY() < max_range[1]) {
			if (p.getX() > max_range[0]) {
				return p.getX() - max_range[0];
			}
			else {
				return min_range[0] - p.getX();
			}
		}
		double dist = 0;
		double anchor[2];
		if (p.getX() > max_range[0]) {
			anchor[0] = max_range[0];
		}
		else {
			anchor[0] = min_range[0];
		}
		if (p.getY() > max_range[1]) {
			anchor[1] = max_range[1];
		}
		else {
			anchor[1] = min_range[1];
		}
		dist = sqrt((anchor[0] - p.getX()) * (anchor[0] - p.getX()) + (anchor[1] - p.getY()) * (anchor[1] - p.getY()));
		return dist;
	}
}
