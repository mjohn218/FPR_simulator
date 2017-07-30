#include "reactions.h"



void perform_association(int p1, Fullmol *bases, Complex *ind_com, Parms &plist, int **crosspart, int **cross_rxn, int **crossint, int **Rlist, int *i_home, double **probvec, int *Ncoup, int **mycoupled, int *p_home, int ci1, int ci2, int it, int *movestat, int *ncross)
{
  
  
  
  int p2 = crosspart[p1][ci1];
  int rxn1 = cross_rxn[p1][ci1];
  
  /*test for being in same complex to avoid moving binding interfaces too far*/
  //flag2=same_complex_test(p1, p2, bases, crossint, i_home, bindrad, rxn1, ci1, ci2);
  cout << "Associate between proteins: p1 " << p1 << ' ' << p2 << " interfaces: " << crossint[p1][ci1] << ' ' << crossint[p2][ci2] << " reaction: " << rxn1 << " pact; " << probvec[p1][ci1] <<  " iter: " << it << endl;
  //				  	assemblyfile<<1<<'\t'<<it*deltat<<'\t'<<p1<<'\t'<<p2<<endl;
  //					cout << bases[p1].xcom << ' ' << bases[p1].ycom << ' ' << bases[p1].zcom <<endl;
  //					cout << bases[p2].xcom << ' ' << bases[p2].ycom << ' ' << bases[p2].zcom <<endl;
  /*Associate proteins, move them to contact and update their free and bound lists*/
  
  
  /*for self binding, use freeleg for clathrin type proteins.*/
  if (bases[p1].protype == plist.pclath && bases[p2].protype ==plist.pclath)
    associate_freelegPBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);
  else
    associate_translate_measurePBCCELL(p1, p2, rxn1, crossint[p1][ci1], crossint[p2][ci2], bases, Rlist, i_home, ind_com, plist);
  
  /*This routine will roate them to align
    their COMs before translating, unless a complex is larger than a single protein.
  */
  
  ap2_coupledrxn_add(rxn1, Ncoup, mycoupled, bases, Rlist, i_home, p_home, p1, p2);
  reflect_complex_rad_rotCELL(p1, bases, ind_com, plist.xboxl, plist.yboxl, plist.zboxl);
  /*Remove p1 and p2 from the list of potential reaction partners so they don't try to
    associate again in this turn.
  */
  remove_reaction_all(p1, p2, ncross, crosspart, probvec, ci1, ci2, cross_rxn, crossint);
  
  /*Since these proteins have moved to associate and taken their complex with them,
    don't allow any proteins in their complex to move again.
  */
  set_movestat_zero(p1, bases, ind_com, movestat);
  
}
