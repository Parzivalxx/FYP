#include "Point.h"

double* FYP::Point::getDataAddr() {
	return values_;
}

bool FYP::Point::equals(Point& point) const {
	if (point.getX() == values_[0] && point.getY() == values_[1] && point.getT() == values_[2]) {
		return true;
	}
	return false;
}
