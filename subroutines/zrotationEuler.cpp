#include "vector_rot_calls.h"

void zrotationEuler(double tx, double ty, double *row) {
	/*Rx=[1 0 0; 0 costhet -sinthet ; 0 sinthet costhet]
	 Ry=[costhet 0 sinthet; 0 1 0 ; -sinthet 0 costhet]
	 Rz=[costhet -sinthet 0 ; sinthet costhet 0 ; 0 0 1]*/

	/*Ryx=[ cosy sinxsiny sinycosx ;
	 0        cosx -sinx ; 
	 -siny cosysinx cosycosx]*/

	/*Rzyx=[  coszcosy, coszsinxsiny -sinzcosx, coszsinycosx+sinzsinx;
	 sinzcosy, sinzsinxsiny+coszcosx, sinzsinycosx-coszsinx;
	 -siny, cosysinx, cosycosx]
	 */

	double sx = sin(tx);
	double cx = cos(tx);
	double sy = sin(ty);
	double cy = cos(ty);

	/*to get the rotated the z component, you only need the last row*/
	//M[6]=-sy;
	//M[7]=cy*sx;
	//M[8]=cy*cx;
	row[0] = -sy;
	row[1] = cy * sx;
	row[2] = cy * cx;

}
