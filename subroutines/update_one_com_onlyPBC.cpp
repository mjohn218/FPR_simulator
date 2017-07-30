#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "reactions.h"
#include "vol_help.h"
#include "rand_gsl.h"
#include "Faddeeva.hh"
#include "utility_calls.h"
#include "vector_rot_calls.h"

void update_one_com_onlyPBC(int c1, Complex *ind_com, Fullmol *bases, Parms plist, int flagbndry) {
	/*Update Complex COM of a single proteins*/
	double totalmassx = 0;
	double totalmassy = 0;
	double totalmassz = 0;
	int s1;
	int mp;
	int i, j;
	int xindex, yindex, zindex;
	ind_com[c1].xcom = 0;
	ind_com[c1].ycom = 0;
	ind_com[c1].zcom = 0;

	/*Size of complex does not change*/
	s1 = ind_com[c1].mysize;

	/*If each dimension spans from quadrant 0 to 2, that's beyond half the boxlength
	 Unless it's a really huge complex (Check for that using radR), that means it's wrapped 
	 around the box.
	 */
	int xnzero = 0;
	int xntwo = 0;
	int ynzero = 0;
	int yntwo = 0;
	int znzero = 0;
	int zntwo = 0;
	//cout <<"Update, complex size: "<<s1<<" complex: "<<c1<<endl;
	for (i = 0; i < s1; i++) {
		mp = ind_com[c1].plist[i];
		//cout <<"in update, curr pos protein: "<<mp<<'\t';
		//cout <<bases[mp].xcom<<' '<<bases[mp].ycom<<' '<<bases[mp].zcom<<endl;
		totalmassx += bases[mp].massx;
		totalmassy += bases[mp].massy;
		totalmassz += bases[mp].massz;
		/*Check for wrapping around box*/
		xindex = int(round((plist.xboxl / 2.0 - bases[mp].xcom) * 2.0 / plist.xboxl)); //0 on right quarter, 1 on middle two quarters, and 2 on left quarter.
		if (xindex == 0)
			xnzero++;
		if (xindex == 2)
			xntwo++;
		yindex = int(round((plist.yboxl / 2.0 - bases[mp].ycom) * 2.0 / plist.yboxl)); //0 on right quarter, 1 on middle two quarters, and 2 on left quarter.
		if (yindex == 0)
			ynzero++;
		if (yindex == 2)
			yntwo++;
		zindex = int(round((plist.zboxl / 2.0 - bases[mp].zcom) * 2.0 / plist.zboxl)); //0 on right quarter, 1 on middle two quarters, and 2 on left quarter.
		if (zindex == 0)
			znzero++;
		if (zindex == 2)
			zntwo++;
		//cout <<"in update com, xindex: "<<xindex<<" bases.xcom: "<<bases[mp].xcom<<" xnzero: "<<xnzero<<" xntwo "<<xntwo<<endl;
		//cout <<"in update com, yindex: "<<yindex<<" bases.xcom: "<<bases[mp].ycom<<" xnzero: "<<ynzero<<" xntwo "<<yntwo<<endl;
		//cout <<"in update com, zindex: "<<zindex<<" bases.xcom: "<<bases[mp].zcom<<" xnzero: "<<znzero<<" xntwo "<<zntwo<<endl;
	}
	int left, right;
	if (xnzero > 0 && xntwo > 0) {
		if (flagbndry == 0) {
			/*wrapped around box*/
			for (i = 0; i < s1; i++) {
				mp = ind_com[c1].plist[i];
				if (bases[mp].xcom < 0)
					ind_com[c1].xcom += (bases[mp].xcom + plist.xboxl) * bases[mp].massx;
				else
					ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
			}
		} else {
			/*This complex is so huge, need another check to determine COM
			 determine if there are proteins that straddle the origin.
			 */

			left = 0;
			right = 0;
			for (i = 0; i < s1; i++) {
				mp = ind_com[c1].plist[i];
				if (bases[mp].xcom < bases[mp].radR && bases[mp].xcom > 0)
					right++;
				else if (bases[mp].xcom > -bases[mp].radR && bases[mp].xcom < 0)
					left++;
			}
			if (left > 0 && right > 0) {
				cout << "Boundary flag=1, tested span of complex, NOT wrapped around box! x" << endl;
				/*calculate regular, non-wrapped COM.*/
				for (i = 0; i < s1; i++) {
					mp = ind_com[c1].plist[i];
					ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
				}
			} else {
				cout << "Boundary flag=1, tested span of complex, DID wrapped around box! x" << endl;
				/*calculate wrapped COM.*/
				for (i = 0; i < s1; i++) {
					mp = ind_com[c1].plist[i];
					if (bases[mp].xcom < 0)
						ind_com[c1].xcom += (bases[mp].xcom + plist.xboxl) * bases[mp].massx;
					else
						ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
				}
			}
		} //checking flagbndry
	} else {

		/*not wrapped around*/
		for (i = 0; i < s1; i++) {
			mp = ind_com[c1].plist[i];
			ind_com[c1].xcom += bases[mp].xcom * bases[mp].massx;
		}
	}
	if (ynzero > 0 && yntwo > 0) {

		if (flagbndry == 0) {
			/*wrapped around box*/
			for (i = 0; i < s1; i++) {
				mp = ind_com[c1].plist[i];
				if (bases[mp].ycom < 0)
					ind_com[c1].ycom += (bases[mp].ycom + plist.yboxl) * bases[mp].massy;
				else
					ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
			}
		} else {
			/*This complex is so huge, need another check to determine COM
			 determine if there are proteins that straddle the origin.
			 */
			left = 0;
			right = 0;
			for (i = 0; i < s1; i++) {
				mp = ind_com[c1].plist[i];
				if (bases[mp].ycom < bases[mp].radR && bases[mp].ycom > 0)
					right++;
				else if (bases[mp].ycom > -bases[mp].radR && bases[mp].ycom < 0)
					left++;
			}
			if (left > 0 && right > 0) {
				cout << "Boundary flag=1, tested span of complex, NOT wrapped around box! y" << endl;
				/*calculate regular, non-wrapped COM.*/
				for (i = 0; i < s1; i++) {
					mp = ind_com[c1].plist[i];
					ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
				}
			} else {
				cout << "Boundary flag=1, tested span of complex, DID wrapped around box! y" << endl;
				/*calculate wrapped COM.*/
				for (i = 0; i < s1; i++) {
					mp = ind_com[c1].plist[i];
					if (bases[mp].ycom < 0)
						ind_com[c1].ycom += (bases[mp].ycom + plist.yboxl) * bases[mp].massy;
					else
						ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
				}

			}
		} //checking flagbndry
	} else {

		/*not wrapped around*/
		for (i = 0; i < s1; i++) {
			mp = ind_com[c1].plist[i];
			ind_com[c1].ycom += bases[mp].ycom * bases[mp].massy;
		}
	}

	if (znzero > 0 && zntwo > 0) {
		if (flagbndry == 0) {
			/*wrapped around box*/
			for (i = 0; i < s1; i++) {
				mp = ind_com[c1].plist[i];
				if (bases[mp].zcom < 0)
					ind_com[c1].zcom += (bases[mp].zcom + plist.zboxl) * bases[mp].massz;
				else
					ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;
			}
		} else {
			/*This complex is so huge, need another check to determine COM
			 determine if there are proteins that straddle the origin.
			 */
			left = 0;
			right = 0;
			for (i = 0; i < s1; i++) {
				mp = ind_com[c1].plist[i];
				if (bases[mp].zcom < bases[mp].radR && bases[mp].zcom > 0)
					right++;
				else if (bases[mp].zcom > -bases[mp].radR && bases[mp].zcom < 0)
					left++;
			}
			if (left > 0 && right > 0) {
				cout << "Boundary flag=1, tested span of complex, NOT wrapped around box! z" << endl;
				/*calculate regular, non-wrapped COM.*/
				for (i = 0; i < s1; i++) {
					mp = ind_com[c1].plist[i];
					ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;
				}
			} else {
				cout << "Boundary flag=1, tested span of complex, DID wrapped around box! z" << endl;
				/*calculate wrapped COM.*/
				for (i = 0; i < s1; i++) {
					mp = ind_com[c1].plist[i];
					if (bases[mp].zcom < 0)
						ind_com[c1].zcom += (bases[mp].zcom + plist.zboxl) * bases[mp].massz;
					else
						ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;
				}

			}

		} //checking flagbndry
	} else {
		/*not wrapped around*/
		for (i = 0; i < s1; i++) {
			mp = ind_com[c1].plist[i];
			ind_com[c1].zcom += bases[mp].zcom * bases[mp].massz;
		}
	}

	ind_com[c1].xcom /= totalmassx;
	ind_com[c1].ycom /= totalmassy;
	ind_com[c1].zcom /= totalmassz;
	//  cout <<"current COM: "<<ind_com[c1].xcom<<' '<<ind_com[c1].ycom<<' '<<ind_com[c1].zcom<<endl;

	/*Now put back in the box*/
	ind_com[c1].xcom -= plist.xboxl * round(ind_com[c1].xcom / plist.xboxl);
	ind_com[c1].ycom -= plist.yboxl * round(ind_com[c1].ycom / plist.yboxl);
	ind_com[c1].zcom -= plist.zboxl * round(ind_com[c1].zcom / plist.zboxl);

}
