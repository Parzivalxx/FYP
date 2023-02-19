#pragma once

#include <iostream>
#include <vector>
#include "../common/Trajectory.h"

namespace FYP {
	class TrajStoreNode {
	public:
		double max_range[3];
		double min_range[3];
		std::vector<Trajectory*>* data_;
		std::vector<int>* children_;
		char* content;
		unsigned int seg_num_;
		bool is_leaf_;
		int id_;
	public:
		TrajStoreNode();
		void updateRange(double x, double y, double t);
		double getXRange() const { return max_range[0] - min_range[0]; };
		double getYRange() const { return max_range[1] - min_range[1]; };
		bool isOverlap(const double query_max[3], const double query_min[3]);
		bool isTemporalOverlap(const double temporal_interval[2]);
		double getL2Distance(Point& p);
	};
}
