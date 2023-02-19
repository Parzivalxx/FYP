#pragma once

#include <vector>
#include <list>
#include "Point.h"
#include "Segment.h"

namespace FYP {
	class Trajectory {
	private:
		unsigned int original_id_;
		unsigned int id_;
		std::vector<Point> points_;

	public:
		Trajectory(){};
		Trajectory(unsigned int id, unsigned int point_num);
		unsigned int getPageSize() const;
		void addPoint(const double x, const double y, const double t);
		void addPoint(const Point& point);
		std::vector<Point>* getPoints() { return &points_; }
		bool isOverlap(const unsigned int dim, const double value) const;
		bool isOverlap(double* query_max, double* query_min);
		bool isTemporalOverlap(double temporal_range[2]);
		void crossSplit(const double cut[2], std::vector<std::list<Trajectory*> >&sub_trajectories);
		bool checkOverlap(Point& p1, Point& p2, double* query_max, double* query_min);
		bool split(const unsigned dim, const double cut, std::vector<unsigned int>& split_pos, std::vector<Point>& split_point);
		int compare(const unsigned int dim, const double value) const; //do not know whether the trajectory is overlapped with the cut
		int fastCompare(const unsigned int dim, const double value)const; //We are sure the trajectory is not overlapped with the cut.
		void clearPoints();
		bool isValid() const;
		unsigned int getSegmentNum() const { return points_.size() - 1; };
		void print() const;
		double getDistanceToPoint(Point& p);
		unsigned int serialize(char* content);
		unsigned int deserialize(char* content);
		unsigned int getByteSize() const;
		void fillSegment(std::vector<Segment*>& segments, int pos);
	};
}
