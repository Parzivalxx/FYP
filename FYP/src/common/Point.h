#pragma once

namespace FYP {
	class Point {
	private:
		double values_[3]; //x, y and t

	public:
		Point(){};
		Point(double x, double y, double t) { values_[0] = x; values_[1] = y; values_[2] = t; }
		Point(double coords[]) { values_[0] = coords[0]; values_[1] = coords[1]; values_[2] = coords[2]; }
		double getX() const { return values_[0]; }
		double getY() const { return values_[1]; }
		double getT() const { return values_[2]; }
		double getValueFromDim(int dim) const { return values_[dim]; }
		void setPoint(double x, double y, double t) {
			values_[0] = x;
			values_[1] = y;
			values_[2] = t;
		}
		void setPoint(Point& p) {
			values_[0] = p.getX();
			values_[1] = p.getY();
			values_[2] = p.getT();
		}
		double* getDataAddr();
		bool equals(Point& p) const;
	};
}
