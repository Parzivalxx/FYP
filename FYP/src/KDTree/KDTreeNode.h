#pragma once

#include <iostream>

namespace FYP {
	class KDTreeNode {
	public:
		double max_range[3];
		double min_range[3];
		unsigned int idx_range_max[3];
		unsigned int idx_range_min[3];
		unsigned int children[2];
		unsigned int parent_;
		vector<int>  data_;
		bool is_leaf_;
		unsigned int depth_;
		unsigned int partition_dim_;
		unsigned int partition_idx_;
		unsigned int segment_num_;
		unsigned int segment_increase_;
		unsigned int overall_segment_;
		bool is_splitted_by_idx_;

	public:
		KDTreeNode();
		bool isOverlap(double* query_max, double* query_min);
		bool isTemporalOverlap(double temporal_range[2]);
		double getL2Distance(Point& p);
		void print()const;

	};
}
