#include "vector_rot_calls.h"

void rotationEulerX(double tx, double *M) {
	/*Rx=[1 0 0; 0 costhet -sinthet ; 0 sinthet costhet]
	 Ry=[costhet 0 sinthet; 0 1 0 ; -sinthet 0 costhet]
	 Rz=[costhet -sinthet 0 ; sinthet costhet 0 ; 0 0 1]*/

	double sx = sin(tx);
	double cx = cos(tx);
	M[0] = 1.0;
	M[1] = 0.0;
	M[2] = 0.0;
	M[3] = 0.0;
	M[4] = cx;
	M[5] = -sx;
	M[6] = 0.0;
	M[7] = sx;
	M[8] = cx;

}
