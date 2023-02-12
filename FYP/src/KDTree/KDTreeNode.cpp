#include "KDTreeNode.h"

FYP::KDTreeNode::KDTreeNode() {

}

bool FYP::KDTreeNode::isTemporalOverlap(double temporal_range[2]) {
	if (max_range[2] < temporal_range[0] || min_range[2] > temporal_range[1]) {
		return false;
	}
	return true;
}

bool FYP::KDTreeNode::isOverlap(double* query_max, double* query_min) {
	for (int i = 0; i < 3; i++) {
		if (query_max[i] < min_range[i] || query_min[i] > max_range[i]) {
			return false;
		}
	}
	return true;
}

void FYP::KDTreeNode::print() const {
	std::cout << "[";
	for (int i = 0; i < 3; i++) {
		std::cout << "(" << min_range[i] << ", " << max_range[i] << ") ";
	}
	std::cout << "]" << std::endl;
}

double FYP::KDTreeNode::getL2Distance(Point& p) {
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
