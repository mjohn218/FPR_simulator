#include "vector_rot_calls.h"

void rotationIntrinsicZYX(double alpha, double beta, double gamma, double *M) {

	/*intrinsic rotation around z, then new y, then new x
	 by first alpha, then beta, then gamma. 
	 So R=Rz(alpha)Ry(beta)Rx(gamma)
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
	/*Rybeta*Rxgamm
	 cb sbsg sbcg
	 0 cg sg
	 -sb cbsg cbcg
	 */

	M[0] = ca * cb;
	M[1] = ca * sb * sg - sa * cg;
	M[2] = ca * sb * cg - sa * sg;
	M[3] = sa * cb;
	M[4] = sa * sb * sg + ca * cg;
	M[5] = sa * sb * cg + ca * sg;
	M[6] = -sb;
	M[7] = sg * cb;
	M[8] = cg * cb;

}
