#include <algorithm>
#include <iostream>

//returns pointer to array
double* findSound(double* c1, double* c2, double* c3, double t1, double t2, double t3) {

	const double v = 1481;
	double c1x = c1[0];
	double c1y = c1[1];
	double c2x = c2[0];
	double c2y = c2[1];
	double c3x = c3[0];
	double c3y = c3[1];

	//make sure c1 is closest to sound
	double tmin = std::min(std::min(t1, t2), t3);
	if (tmin == t2) {
		std::swap(t1, t2);
		std::swap(c1x, c2x);
		std::swap(c1y, c2y);
	}
	if (tmin == t3) {
		std::swap(t1, t3);
		std::swap(c1x, c3x);
		std::swap(c1y, c3y);
	}

	double d1 = v*t1;
	double t12 = t2-t1;
	double t13 = t3-t1;
	double A = 2*c1x-2*c2x;
	double B = 2*c1y-2*c2y;
	double C = (v*t12)*(v*t12) + 2*v*t12*d1 + c1x*c1x + c1y*c1y - c2x*c2x - c2y*c2y;
	double D = 2*c1x-2*c3x;
	double E = 2*c1y-2*c3y;
	double F = (v*t13)*(v*t13) + 2*v*t13*d1 + c1x*c1x + c1y*c1y - c3x*c3x - c3y*c3y;

	double soundY = (C*D-A*F)/(B*D-A*E);
	double soundX = (C-B*soundY)/A;

	double* sound = new double[2];
	sound[0] = soundX;
	sound[1] = soundY;
	return sound;
}

int main() {

	double* sound;
	double c1[2] = {.5,0};
	double c2[2] = {-.5,0};
	double c3[2] = {0,1};
	double t1 = .00979068;
	double t2 = .01026799;
	double t3 = .00956097;
	sound = findSound(c1, c2, c3, t1, t2, t3);
	std::cout << sound[0] << ' ' << sound[1] << std::endl;
	return 0;
}
