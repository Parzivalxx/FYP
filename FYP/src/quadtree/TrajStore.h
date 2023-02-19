#pragma once

#include "TrajStoreNode.h"

namespace FYP {

	class TrajStore {
	private:
		unsigned int page_size_ = 1024;
		double query_spatial_domain[2] = { 100.0, 100.0 };
		double epsilon = 0;
		int fanout = 50;
		double stop_increase_ratio = 0;
		std::vector<TrajStoreNode*> tree_nodes_ = {};

		int node_access = 0;
		int leaf_access = 0;
		int block_access = 0;

	public:
		~TrajStore();
		void setTSStopIncreaseRatio(double ratio);
		void setTSPageSize(unsigned int page_size);
		void setTSQueryPara(double qx, double qy);
		void setTSEpsilon(double epsilon);
		void setTSFanout(int fanout);
		int getTSNodeAccess();
		int getTSLeafAccess();
		int getTSBlockAccess();
		void initialize(const std::string traj_file_name);
		bool stopCondition(int node_id);
		void finalizeNode(int node_id);
		void rangeQuery(double query_max[3], double query_min[3]);
		void kNNQuery(int k, double temporal_interval[2], double query_location[2]);
	};
}
