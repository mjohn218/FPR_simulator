#include "vector_rot_calls.h"

void BodyFixedRot() {
	/*To rotate around the axes of the body-fixed axes in the molecule,
	 need rotation matrixes around local z, y, and x axes.
	 One can then perform a series of intrinsic rotations using these matrices,
	 R=Z(sqrt(Dz 2dt)) Y( sqrt(Dy 2dt) X( sqrt(Dx 2dt))

	 */

}
