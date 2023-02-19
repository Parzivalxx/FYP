#pragma once

#include "../quadtree/TrajStore.h"

namespace FYP {
	double getPDistance(double x, double y, double x1, double y1, double x2, double y2);
	double getCrossProduct(double* v1, double* v2);
	bool isLineRectangleOverlap(double x_start, double y_start, double x_end, double y_end, double left, double bottom, double right, double top);
}
