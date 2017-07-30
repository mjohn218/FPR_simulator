#include "reactions.h"
#include "utility_calls.h"

void update_rates(Protein *wholep, Parms &plist, int **Rlist, double *bindrad, double *kr, int *rxtype, double *Kd, int *phome) {

	/*Use this if you are reading in from the reaction file the values of the experimental on rates, rather that the intrinsic
	 activation rates, because the intrinsic activation rate is used in the GF probability of reaction occuring.
	 The different between the experimental on rate, k_on, and intrinsic activation rate, k_a, is defined by the equation
	 1/k_on= 1/k_a + 1/k_D. k_on is the experimental on rate that accounts for both
	 the intrinsic reaction rate k_a, and the diffusion to contact, k_D=4pi*sigma*D_tot
	 
	 Same must be done for k_off, where 1/k_off=1/k_d + Keq/k_D, and k_d is the intrinsic off rate and Keq is the equilibration
	 constant.
	 */
	int p1, p2;
	int j;
	double fourpi = 4.0 * M_PI;
	double kdiff, Dtot;
	double kactinv;
	int i1, i2;
	for (j = 0; j < plist.Nrxn; j++) {
		if (rxtype[j] == 0) {
			/*Then this is a binding reaction*/
			i1 = Rlist[j][0];
			i2 = Rlist[j][1];
			p1 = phome[i1];
			p2 = phome[i2];
			Dtot = wholep[p1].Dx + wholep[p2].Dx;
			kdiff = fourpi * Dtot * bindrad[j];
			kactinv = 1.0 / kr[j] - 1.0 / (kdiff);
			if (kr[j] != 0)
				kr[j] = 1.0 / kactinv;
			cout << "reaction: " << j << " new ka: " << kr[j] << endl;

		} else if (rxtype[j] == 1) {
			i1 = Rlist[j][1];
			i2 = Rlist[j][2];
			p1 = phome[i1];
			p2 = phome[i2];

			Dtot = wholep[p1].Dx + wholep[p2].Dx;
			kdiff = fourpi * Dtot * bindrad[j]; //nm^3/us
			Kd[j] = Kd[j] * 1E-6 / 1E24 * 6.022E23; //to put in units of molecule/nm^3
			kactinv = 1.0 / kr[j] - 1.0 / (Kd[j] * 1E6 * kdiff); //to put in units of s^-1 
			if (kr[j] != 0)
				kr[j] = 1.0 / kactinv;

			cout << "unbind reaction: " << j << " new kd: " << kr[j] << endl;

		}
		/*zeroth order doesn't change*/

	}

}
