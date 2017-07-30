#include "2Drelated.h"

double fpir_cum(double x, void *p) {

	double alp, bet, tet, P, T, h, f;
	f_params &params = *reinterpret_cast<f_params *>(p);

	h = (2.0 * M_PI * params.a * params.D);

	if (params.k < 1E23) {
		alp = h * x * j1(x * params.a) + params.k * j0(x * params.a);
		bet = h * x * y1(x * params.a) + params.k * y0(x * params.a);
		tet = sqrt(alp * alp + bet * bet);

		P = (j1(x * params.r) *params.r* bet - y1(x * params.r) *params.r* alp) / tet;
		T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;

		f =  exp(-params.D * params.t * x * x) * P * T / (2.0 * M_PI);
	} else { //absorbing boundary conditions
		alp = j0(x * params.a);
		bet = y0(x * params.a);
		tet = sqrt(alp * alp + bet * bet);

		P = (j0(x * params.r) * bet - y0(x * params.r) * alp) / tet;
		T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;

		f = x * exp(-params.D * params.t * x * x) * P * T / (2.0 * M_PI);
	}

	return f;
}
