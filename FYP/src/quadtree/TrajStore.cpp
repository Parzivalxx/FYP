#include "TrajStore.h"
#include "../common/utils.h"

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <queue>
#include <set>


FYP::TrajStore::~TrajStore() {
	for (TrajStoreNode* node : tree_nodes_) {
		delete node;
	}
	tree_nodes_.clear();
}

void FYP::TrajStore::setTSStopIncreaseRatio(double ratio) {
	stop_increase_ratio = ratio;
}

void FYP::TrajStore::setTSPageSize(unsigned int page_size) {
	page_size_ = page_size;
}

void FYP::TrajStore::setTSQueryPara(double qx, double qy) {
	query_spatial_domain[0] = qx;
	query_spatial_domain[1] = qy;
}

void FYP::TrajStore::setTSEpsilon(double epsilon) {
	epsilon = epsilon;
}

void FYP::TrajStore::setTSFanout(int fanout) {
	fanout = fanout;
}

int FYP::TrajStore::getTSNodeAccess() {
	return node_access;
}

int FYP::TrajStore::getTSLeafAccess() {
	return leaf_access;
}

int FYP::TrajStore:: getTSBlockAccess() {
	return block_access;
}

void FYP::TrajStore::initialize(const std::string traj_file_name) {
	std::ifstream ifs(traj_file_name, std::ios::in);
	unsigned int total_traj = 0;
	unsigned int point_num = 0;
	ifs >> total_traj;
	std::cout << total_traj << " trajectories" << std::endl;
	TrajStoreNode* root = new TrajStoreNode();
	root->data_ = new std::vector<Trajectory*>(total_traj, nullptr);
	root->id_ = tree_nodes_.size();
	std::cout << "root node constructed" << std::endl;
	for (unsigned int i = 0; i < total_traj; i++) {
		ifs >> point_num;
		//cout << i << " " << point_num << endl;
		root->data_->at(i) = new Trajectory(i, point_num);
		for (unsigned int point_id = 0; point_id < point_num; point_id += 1) {
			double x, y, t;
			ifs >> x >> y >> t;
			root->updateRange(x, y, t);
			root->data_->at(i)->getPoints()->at(point_id).setPoint(x, y, t);
		}
	}
	ifs.close();
	std::cout << "trajectories loaded" << std::endl;
	std::getchar();
	tree_nodes_.push_back(root);
	std::list<TrajStoreNode*> queue;
	queue.push_back(root);
	while (!queue.empty()) {
		TrajStoreNode* node = queue.front();
		queue.pop_front();
		bool success = stopCondition(node->id_);

		if (success) {
			for (int child_id : *(node->children_)) {
				TrajStoreNode* child = tree_nodes_[child_id];
				queue.push_back(child);
			}
		}
		else {
			finalizeNode(node->id_);
		}
		//cout <<"success: "<< success << " queue.size(): " << queue.size() << endl;
	}
}

void FYP::TrajStore::finalizeNode(int node_id) {
	TrajStoreNode* node = tree_nodes_[node_id];
	if (node->data_ == nullptr) {
		std::cout << "Node " << node_id << " does not have any data yet." << std::endl;
		exit(0);
	}
	int total_segments_num = 0;
	for (Trajectory* traj : *(node->data_)) {
		total_segments_num += traj->getPoints()->size() - 1;
	}
	node->is_leaf_ = false;
	std::vector<Segment*> all_segments(total_segments_num, nullptr);
	int start_pos = 0;
	for (Trajectory* traj : *(node->data_)) {
		traj->fillSegment(all_segments, start_pos);
		start_pos += traj->getPoints()->size() - 1;
	}
	sort(all_segments.begin(), all_segments.end(), [](Segment* lhs, Segment* rhs) {return lhs->getEnd()->getT() < rhs->getEnd()->getT(); });
	unsigned int size = all_segments.size() * (sizeof(unsigned int) + sizeof(double) * 3 * 2);
	int child_num = (int)ceil(1.0 * all_segments.size() / fanout);
	node->children_ = new std::vector<int>(child_num);
	for (int c_id = 0; c_id < child_num; c_id++) {

		TrajStoreNode* child = new TrajStoreNode();
		child->id_ = tree_nodes_.size();
		child->is_leaf_ = true;
		node->children_->at(c_id) = child->id_;
		child->seg_num_ = all_segments.size() - c_id * fanout;
		if (child->seg_num_ > fanout) {
			child->seg_num_ = fanout;
		}
		child->content = new char[child->seg_num_ * (sizeof(unsigned int) + sizeof(double) * 6)];
		unsigned int start_idx = 0;
		for (int i = 0; i < child->seg_num_; i++) {
			Segment* seg = all_segments[c_id * fanout + i];
			child->updateRange(seg->getStart()->getX(), seg->getStart()->getY(), seg->getStart()->getT());
			child->updateRange(seg->getEnd()->getX(), seg->getEnd()->getY(), seg->getEnd()->getT());
			int len = seg->serialize(child->content + start_idx);
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



bool FYP::TrajStore::stopCondition(int node_id) {
	TrajStoreNode* node = tree_nodes_[node_id];
	//cout << "node range: " << node->min_range[0] << " " << node->max_range[0] << ", " << node->min_range[1] << " " << node->max_range[1] << endl;
	unsigned int byte_size = 0;
	unsigned int seg_num_before_split = 0;
	unsigned int seg_num_after_split = 0;
	//cout << "trajectory: " << endl;
	for (int i = 0; i < node->data_->size(); i++) {
		byte_size += node->data_->at(i)->getByteSize();
		//node->data_->at(i)->Print();
		seg_num_before_split += node->data_->at(i)->getSegmentNum();
	}
	if (byte_size <= page_size_) {
		//no need to split
		return false;
	}
	double cost_not_split = (query_spatial_domain[0] + node->getXRange()) * (query_spatial_domain[1] + node->getYRange()) * ceil(1.0 * byte_size / page_size_);

	TrajStoreNode* children[4];
	for (int i = 0; i < 4; i++) {
		children[i] = new TrajStoreNode();
		children[i]->data_ = new std::vector<Trajectory*>();
		children[i]->data_->reserve(node->data_->size());
	}
	std::vector<std::list<Trajectory*> > sub_trajectories(4);
	double cross_cuts[2] = { 0.5 * (node->max_range[0] + node->min_range[0]), 0.5 * (node->max_range[1] + node->min_range[1]) };
	for (int traj_id = 0; traj_id < node->data_->size(); traj_id++) {
		Trajectory* traj = node->data_->at(traj_id);
		traj->crossSplit(cross_cuts, sub_trajectories);
		for (int child_id = 0; child_id < 4; child_id++) {
			for (Trajectory* subtraj : sub_trajectories[child_id]) {
				//subtraj->Print();
				children[child_id]->data_->push_back(subtraj);
				seg_num_after_split += subtraj->getSegmentNum();
			}
			sub_trajectories[child_id].clear();
		}
	}
	double cost_after_split = 0;
	unsigned int split_byte_sizes[4] = { 0, 0, 0, 0 };
	for (int child_id = 0; child_id < 4; child_id++) {
		for (Trajectory* traj : *(children[child_id]->data_)) {
			split_byte_sizes[child_id] += traj->getByteSize();
			for (Point& p : *(traj->getPoints())) {
				children[child_id]->updateRange(p.getX(), p.getY(), p.getT());
			}
		}
		cost_after_split += (query_spatial_domain[0] + 0.5 * node->getXRange()) * (query_spatial_domain[1] + 0.5 * node->getYRange()) * ceil(1.0 * split_byte_sizes[child_id] / page_size_);
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
		node->children_ = new std::vector<int>();
		for (int i = 0; i < 4; i++) {
			if (children[i]->data_->empty()) {

				if (children[i]->data_ == nullptr) {
					std::cout << "null" << std::endl;
				}
				if (children[i] == nullptr) {
					std::cout << "null2" << std::endl;
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

void FYP::TrajStore::kNNQuery(int k, double temporal_interval[2], double coord[2]) {
	node_access = 0;
	leaf_access = 0;
	block_access = 0;
	Point q;
	q.setPoint(coord[0], coord[1], 0.0);
	TrajStoreNode* root = tree_nodes_[0];
	std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int> >, std::greater<std::pair<double, int> > > pqueue;
	if (root->isTemporalOverlap(temporal_interval)) {
		pqueue.emplace(root->getL2Distance(q), 0);
	}
	std::vector<std::pair<double, int> > results;
	while (!pqueue.empty()) {
		if (results.size() == k && results.back().first < pqueue.top().first) {
			break;
		}
		std::pair<double, int> top = pqueue.top();
		pqueue.pop();
		TrajStoreNode* node = tree_nodes_[top.second];
		if (node->is_leaf_) {
			leaf_access += 1;
			node_access += 1;
			block_access += (int)ceil(1.0 * node->seg_num_ * (sizeof(unsigned int) + sizeof(double) * 6) / page_size_);
			Segment seg;
			unsigned int pos = 0;
			for (int i = 0; i < node->seg_num_; i++) {
				unsigned int l = seg.deserialize(node->content + pos);
				pos += l;
				if (!seg.isTemporalOverlap(temporal_interval)) {
					continue;
				}
				double distance = getPDistance(q.getX(), q.getY(), seg.getStart()->getX(), seg.getStart()->getY(), seg.getEnd()->getX(), seg.getEnd()->getY());
				int original_id = seg.getOriginalId();
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
				std::cout << "node has not children" << std::endl;
				exit(0);
			}
			for (int c_id : *(node->children_)) {
				TrajStoreNode* child = tree_nodes_[c_id];
				if (child->isTemporalOverlap(temporal_interval)) {
					double d = child->getL2Distance(q);
					pqueue.emplace(d, c_id);
				}
			}
		}
	}
}

void FYP::TrajStore::rangeQuery(double query_max[3], double query_min[3]) {
	node_access = 0;
	leaf_access = 0;
	block_access = 0;
	std::list<int> node_queue;
	TrajStoreNode* root = tree_nodes_[0];
	if (root->isOverlap(query_max, query_min)) {
		node_queue.push_back(0);
	}
	std::set<int> results;
	while (!node_queue.empty()) {
		int node_id = node_queue.front();
		node_queue.pop_front();
		TrajStoreNode* node = tree_nodes_[node_id];
		if (node->is_leaf_) {
			Segment seg;
			unsigned int pos = 0;
			for (int seg_id = 0; seg_id < node->seg_num_; seg_id++) {
				unsigned int l = seg.deserialize(node->content + pos);
				pos += l;
				if (seg.isOverlap(query_max, query_min)) {
					results.insert(seg.getOriginalId());
				}
			}
			node_access += 1;
			leaf_access += 1;
			block_access += (int)ceil(1.0 * node->seg_num_ * (sizeof(unsigned int) + sizeof(double) * 6) / page_size_);
		}
		else {
			if (node->children_ == nullptr) {
				std::cout << "node has no children" << std::endl;
				exit(0);
			}
			for (int c_id : *(node->children_)) {
				TrajStoreNode* child = tree_nodes_[c_id];
				if (child->isOverlap(query_max, query_min)) {
					node_queue.push_back(c_id);
				}
			}
			node_access += 1;
		}
	}
	std::cout << results.size() << " trajectories retrieved" << std::endl;
	/*for (int i : results) {
		cout << i << " ";
	}
	cout << endl;*/
}
