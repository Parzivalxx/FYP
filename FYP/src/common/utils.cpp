#include "utils.h"
#include <math.h>

double FYP::getPDistance(double x, double y, double x1, double y1, double x2, double y2) {
	double a = x - x1;
	double b = y - y1;
	double c = x2 - x1;
	double d = y2 - y1;

	double dot = a * c + b * d;
	double len_sq = c * c + d * d;
	double param = -1;
	if (len_sq != 0) {
		param = dot / len_sq;
	}
	double xx, yy;
	if (param < 0) {
		xx = x1;
		yy = y1;
	}
	else if (param > 1) {
		xx = x2;
		yy = y2;
	}
	else {
		xx = x1 + param * c;
		yy = y1 + param * d;
	}
	double dx = x - xx;
	double dy = y - yy;
	return sqrt(dx * dx + dy * dy);
}

double FYP::getCrossProduct(double* v1, double* v2) {
	return v1[0] * v2[1] - v2[0] * v1[1];
}

bool FYP::isLineRectangleOverlap(double x_start, double y_start, double x_end, double y_end, double left, double bottom, double right, double top) {
	//cout << x_start << " " << y_start << " " << x_end << " " << y_end << " " << top << " " << bottom << " " << left << " " << right << endl;
	if (left <= x_start && x_start <= right && bottom <= y_start && y_start <= top) {
		//the start point is inside the rectangle;
		return true;
	}
	if (left <= x_end && x_end <= right && bottom <= y_end && y_end <= top) {
		//the end point is inside the rectangle;
		return true;
	}
	double v1[2]{};
	double v2[2]{};
	double v3[2]{};
	double v4[2]{};
	double v5[2]{};
	double v6[2]{};
	v1[0] = x_start - right;
	v1[1] = y_start - bottom;
	v2[0] = left - right;
	v2[1] = top - bottom;
	v3[0] = x_end - right;
	v3[1] = y_end - bottom;

	v4[0] = left - x_end;
	v4[1] = top - y_end;
	v5[0] = x_start - x_end;
	v5[1] = y_start - y_end;
	v6[0] = right - x_end;
	v6[1] = bottom - y_end;

	/*cout << left << " " << right << " " << bottom << " " << top << endl;
	cout << x_start << " " << y_start << " " << x_end << " " << y_end << endl;
	cout << "cross product: " << CrossProduct(v1, v2) << " " << CrossProduct(v3, v2) << endl;
	*/
	if ((getCrossProduct(v1, v2) * getCrossProduct(v3, v2) <= 0) && (getCrossProduct(v4, v5) * getCrossProduct(v6, v5) <= 0)) {
	//segment intersect with the diagonal line1
		return true;
	}
	v1[0] = x_start - left;
	v1[1] = y_start - bottom;
	v2[0] = right - left;
	v2[1] = top - bottom;
	v3[0] = x_end - left;
	v3[1] = y_end - bottom;
	v4[0] = left - x_end;
	v4[1] = bottom - y_end;
	v5[0] = x_start - x_end;
	v5[1] = y_start - y_end;
	v6[0] = right - x_end;
	v6[1] = top - y_end;
	//cout << "cross product: " << CrossProduct(v1, v2) << " " << CrossProduct(v3, v2) << endl;
	if ((getCrossProduct(v1, v2) * getCrossProduct(v3, v2) <= 0) && (getCrossProduct(v4, v5) * getCrossProduct(v6, v5) <= 0)) {
		//segment intersect with the diagonal line2
		return true;
	}
	return false;
}
