#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>

using namespace std;

double Vcap_three_dim(int flagx, int flagy, int flagz, double dx, double dy, double dz, double R);
double vcap_overlap(double dx, double dy, double R);
double Vcap_one_dim(int flagx, int flagy, double dx, double dy, double dz, double R);
double vcap_three(double dx, double dy, double dz, double R);
double Vcap_two_dim(int flagx, int flagy, double dx, double dy, double dz, double R);
double integrate_atanx2(double p, double R, double q1, double q2);
double integrate_atan(double p, double R, double q1, double q2);
double integrate_atanx2inv(double p, double R, double q1, double q2);
double integrate_ataninv(double p, double R, double q1, double q2);
double sphere_cap(double R, double h);
double area_slice(double dx, double dy, double dz, double R);
void runmctest(double r1, double xcom, double ycom, double zcom, double Vcalc, double box_x);
void numerical_test(double dx, double dy, double dz, double r1);

