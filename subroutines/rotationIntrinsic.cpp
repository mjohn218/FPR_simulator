#include "vector_rot_calls.h"

void rotationIntrinsic(double alpha, double beta, double gamma, double *M) {
	/*Don't use this for rotational diffusion!
	 
	 */
	/*intrinsic rotation around z, then new x, then new z
	 by first alpha, then beta, then gamma. 
	 So R=Rz(alpha)Rx(beta)Rz(gamma)
	 reverse order of an extrinsic rotation
	 
	 */
	/*THESE MATRIX DEFINITIONS REQUIRE THE LOCAL FRAME INITIALLY
	 ALIGNS WITH THE x,y.z AXES!*/
	/*Rx=[1 0 0; 
	 0 costhet -sinthet ; 
	 0 sinthet costhet]
	 Ry=[costhet 0 sinthet; 0 1 0 ; -sinthet 0 costhet]
	 Rz=[costhet -sinthet 0 ; sinthet costhet 0 ; 0 0 1]*/
	double sa = sin(alpha);
	double ca = cos(alpha);
	double sb = sin(beta);
	double cb = cos(beta);
	double sg = sin(gamma);
	double cg = cos(gamma);
	/*Rxbeta*Rzgamm
	 [cg -sg 0
	 cb*sg cb*cg -sb
	 sb*sg sb*cg cb
	 */
	/*Rzalpha*Rxbeta*Rzgamma
	 ca*cg-sa*cb*sg  -ca*sg-sa*cb*cg  sa*sb
	 sa*cg+ca*cb*sg   -sa*sg+ca*cb*cg -ca*sb
	 sb*sg   sb*cg   cb 
	 */
	M[0] = ca * cg - sa * cb * sg;
	M[1] = -ca * sg - sa * cb * cg;
	M[2] = sa * sb;
	M[3] = sa * cg + ca * cb * sg;
	M[4] = -sa * sg + ca * cb * cg;
	M[5] = -ca * sb;
	M[6] = sb * sg;
	M[7] = sb * cg;
	M[8] = cb;

}
