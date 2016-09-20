#include <algorithm>
#include <iostream>
#include <math.h>

//hyperbola function
double hbFunc(double* p, double* c0, double* c1, double d) {

	double x = p[0];
	double y = p[1];
	double c0x = c0[0];
	double c0y = c0[1];
	double c1x = c1[0];
	double c1y = c1[1];

	return sqrt((x-c1x)*(x-c1x)+(y-c1y)*(y-c1y))-sqrt((x-c0x)*(x-c0x)+(y-c0y)*(y-c0y))-d;
}

//hyperbola derivative x
double hbDerivX(double* p, double* c0, double* c1, double d) {

	double x = p[0];
	double y = p[1];
	double c0x = c0[0];
	double c0y = c0[1];
	double c1x = c1[0];
	double c1y = c1[1];

	return (x-c1x)/sqrt((x-c1x)*(x-c1x)+(y-c1y)*(y-c1y))-(x-c0x)/sqrt((x-c0x)*(x-c0x)+(y-c0y)*(y-c0y));
}

//hyperbola derivative y
double hbDerivY(double* p, double* c0, double* c1, double d) {

	double x = p[0];
	double y = p[1];
	double c0x = c0[0];
	double c0y = c0[1];
	double c1x = c1[0];
	double c1y = c1[1];

	return (y-c1y)/sqrt((x-c1x)*(x-c1x)+(y-c1y)*(y-c1y))-(y-c0y)/sqrt((x-c0x)*(x-c0x)+(y-c0y)*(y-c0y));
}

//intersection of two lines
double* findIntersect(double s1, double s2, double p1x, double p1y, double p2x, double p2y) {

	double* intersect = new double[2];
	intersect[0] = (s1*p1x-p1y-s2*p2x+p2y)/(s1-s2);
	intersect[1] = s1*(intersect[0]-p1x)+p1y;
	return intersect;
}

//main function
double* findSound(double* c0, double* c1, double* c2, double t0, double t1, double t2) {

	const double v = 1481;
	double c0x = c0[0];
	double c0y = c0[1];
	double c1x = c1[0];
	double c1y = c1[1];
	double c2x = c2[0];
	double c2y = c2[1];

	//make sure c0 is closest to sound
	double tmin = std::min(std::min(t0, t1), t2);
	if (tmin == t1) {
		std::swap(t0, t1);
		std::swap(c0x, c1x);
		std::swap(c0y, c1y);
		std::swap(c0, c1);
	}
	if (tmin == t2) {
		std::swap(t0, t2);
		std::swap(c0x, c2x);
		std::swap(c0y, c2y);
		std::swap(c0, c2);
	}

	//find asymptotes of hyperbolas for first estimate
	//h stands for hyperbola value
	//lots of math and repetition here
	double c1h = sqrt((c0x-c1x)*(c0x-c1x)+(c0y-c1y)*(c0y-c1y))/2;
	double a1h = v*(t1-t0)/2;
	double d1 = 2*a1h;
	double b1h_1 = sqrt(c1h*c1h-a1h*a1h);
	double b1h_2 = -b1h_1;
	double theta1 = atan((c1y-c0y)/(c1x-c0x));
	double slope1_1 = (a1h*sin(theta1)+b1h_1*cos(theta1))/(a1h*cos(theta1)-b1h_1*sin(theta1));
	double slope1_2 = (a1h*sin(theta1)+b1h_2*cos(theta1))/(a1h*cos(theta1)-b1h_2*sin(theta1));
	double p1x = (c0x+c1x)/2;
	double p1y = (c0y+c1y)/2;
	//std::cout << "c="<<c1h<<"\na="<<a1h<<"\nb="<<b1h_1<<"\ntheta="<<theta1<<"\nslope1="<<slope1_1<<"\nslope2="<<slope1_2<<"\npoint=("<<p1x<<","<<p1y<<")\n";

	double c2h = sqrt((c0x-c2x)*(c0x-c2x)+(c0y-c2y)*(c0y-c2y))/2;
	double a2h = v*(t2-t0)/2;
	double d2 = 2*a2h;
	double b2h_1 = sqrt(c2h*c2h-a2h*a2h);
	double b2h_2 = -b2h_1;
	double theta2 = atan((c2y-c0y)/(c2x-c0x));
	double slope2_1 = (a2h*sin(theta2)+b2h_1*cos(theta2))/(a2h*cos(theta2)-b2h_1*sin(theta2));
	double slope2_2 = (a2h*sin(theta2)+b2h_2*cos(theta2))/(a2h*cos(theta2)-b2h_2*sin(theta2));
	double p2x = (c0x+c2x)/2;
	double p2y = (c0y+c2y)/2;
	//std::cout << slope1_1 << '\n';
	//std::cout << slope1_2 << '\n';
	//std::cout << slope2_1 << '\n';
	//std::cout << slope2_2 << '\n';

	double* intersect1 = findIntersect(slope1_1, slope2_1, p1x, p1y, p2x, p2y);
	double* intersect2 = findIntersect(slope1_2, slope2_1, p1x, p1y, p2x, p2y);
	double* intersect3 = findIntersect(slope1_1, slope2_2, p1x, p1y, p2x, p2y);
	double* intersect4 = findIntersect(slope1_2, slope2_2, p1x, p1y, p2x, p2y);
	//std::cout << "first intersection is: ("<<intersect1[0]<<','<<intersect1[1]<<")\n";
	//std::cout << "second intersection is: ("<<intersect2[0]<<','<<intersect2[1]<<")\n";
	//std::cout << "third intersection is: ("<<intersect3[0]<<','<<intersect3[1]<<")\n";
	//std::cout << "fourth intersection is: ("<<intersect4[0]<<','<<intersect4[1]<<")\n";

	double val1 = hbFunc(intersect1, c0, c1, d1)*hbFunc(intersect1, c0, c2, d2);
	double val2 = hbFunc(intersect2, c0, c1, d1)*hbFunc(intersect2, c0, c2, d2);
	double val3 = hbFunc(intersect3, c0, c1, d1)*hbFunc(intersect3, c0, c2, d2);
	double val4 = hbFunc(intersect4, c0, c1, d1)*hbFunc(intersect4, c0, c2, d2);
	//std::cout << "val1="<<val1<<'\n';
	//std::cout << "val2="<<val2<<'\n';
	//std::cout << "val3="<<val3<<'\n';
	//std::cout << "val4="<<val4<<'\n';

	double* guess;
	double guessVal = std::min(std::min(std::min(val1,val2),val3),val4);
	if (guessVal == val1) guess = intersect1;
	if (guessVal == val2) guess = intersect2;
	if (guessVal == val3) guess = intersect3;
	if (guessVal == val4) guess = intersect4;

	//newton's method to improve guess
	double guessDiff = 1000;
	double nextGuess[2];
	double A,B,C,D,E,F;
	while (guessDiff > .00000001) {
		//system of equations:
		//Ax+By=C
		//Dx+Ey=F
		//std::cout << guess[0] << ' ' << guess[1] << '\n';
		A = hbDerivX(guess, c0, c1, d1);
		B = hbDerivY(guess, c0, c1, d1);
		C = guess[0]*hbDerivX(guess, c0, c1, d1)+guess[1]*hbDerivY(guess, c0, c1, d1)-hbFunc(guess, c0, c1, d1);
		D = hbDerivX(guess, c0, c2, d2);
		E = hbDerivY(guess, c0, c2, d2);
		F = guess[0]*hbDerivX(guess, c0, c2, d2)+guess[1]*hbDerivY(guess, c0, c2, d2)-hbFunc(guess, c0, c2, d2);

		nextGuess[0] = (C*E-B*F)/(A*E-B*D);
		nextGuess[1] = (C*D-A*F)/(B*D-A*E);
		//guessDiff = sqrt((guess[0]-nextGuess[0])*(guess[0]-nextGuess[0])+(guess[1]-nextGuess[1])*(guess[0]-nextGuess[0]));
		guessDiff = fabs(guess[0]-nextGuess[0])+fabs(guess[1]-nextGuess[1]);
		guess[0] = nextGuess[0];
		guess[1] = nextGuess[1];
		//std::cout << "guessDiff="<<guessDiff<<'\n';
	}

	double* sound = new double[2];
	sound[0] = guess[0];
	sound[1] = guess[1];
	return sound;

}

int main() {

	double* sound;
	double c0[2] = {.5,0};
	double c1[2] = {-.5,0};
	double c2[2] = {0,1};
	double t0 = 1.00070702;
	double t1 = 1.00026347;
	double t2 = 1.0;

//other random test cases
/*
	double t0 = .00979068;
	double t1 = .01026799;
	double t2 = .00956097;
*/
/*
	t0 = 3.603774592;
	t1 = 3.604121052;
	t2 = 3.604529510;
*/
	
	sound = findSound(c0, c1, c2, t0, t1, t2);
	std::cout << sound[0] << ' ' << sound[1] << std::endl;
	return 0;
}
