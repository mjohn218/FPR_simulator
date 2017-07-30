#include "vector_rot_calls.h"

void rotationEulerZXZ(double alpha, double beta, double gamma, double *M) {
	/*Extrinsic around z by alpha, then x by beta. then z by gamma
	 */
	double sa = sin(alpha);
	double ca = cos(alpha);
	double sb = sin(beta);
	double cb = cos(beta);
	double sg = sin(gamma);
	double cg = cos(gamma);

	/*Rzgamma*Rxbeta*Rzalpha
	 
	 */
	M[0] = cg * ca - sg * cb * sa;
	M[1] = -cg * sa - sg * cb * ca;
	M[2] = sg * sb;
	M[3] = sg * ca + cg * cb * sa;
	M[4] = -sg * sa + cg * cb * ca;
	M[5] = -cg * sb;
	M[6] = sb * sa;
	M[7] = sb * ca;
	M[8] = cb;

}
