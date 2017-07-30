#include "reactions.h"

void remove_reaction_all_ncom(int p1, int p2, int *ncross, int **crosspart, double **prob, int **cross_rxn, int **cross_int, int *ncrosscom, Fullmol *bases)
{ 
  
  /*remove proteins p1 and p2 from the list of their potential reaction partners*/
  /*After doing this, you won't avoid overlap with any of your current close enough neighbors, and they won't avoid overlap with you
    having been deleted from their list.
   */
  int pskip;
  int i, j;
  double pmatch;
  
  for(i=0;i<ncross[p1];i++){
    
    pskip=crosspart[p1][i];
    if(pskip!=p2){
      /*this criteria is needed, otherwise you might delete p1 from p2's list if they interact in 2 ways*/
      
      /*Need to delete protein p1 from the pskip's list    */
      
      for(j=0;j<ncross[pskip];j++){
	if(crosspart[pskip][j]==p1){
	  crosspart[pskip][j]=crosspart[pskip][ncross[pskip]-1];
	  prob[pskip][j]=prob[pskip][ncross[pskip]-1];
	  cross_rxn[pskip][j]=cross_rxn[pskip][ncross[pskip]-1];
	  cross_int[pskip][j]=cross_int[pskip][ncross[pskip]-1];
	  ncross[pskip]-=1;
	  ncrosscom[bases[pskip].mycomplex]-=1;
	  j-=1;
	  //pskip and p1 may interact in more than one way, delete all times
	}
      }
      
      
    }

  }
  
  int p3;
  /*p2 has reacted, so if it has multiple partners, remove it from their lists,
    same as above.
  */

  if(ncross[p2]>1){
    
    for(i=0;i<ncross[p2];i++){
      p3=crosspart[p2][i];
      if(p3!=p1){
	/*don't delete p2 from p1's list of partners*/
	
	for(j=0;j<ncross[p3];j++){
	  if(crosspart[p3][j]==p2){
	    crosspart[p3][j]=crosspart[p3][ncross[p3]-1];
	    prob[p3][j]=prob[p3][ncross[p3]-1];
	    cross_rxn[p3][j]=cross_rxn[p3][ncross[p3]-1];
	    cross_int[p3][j]=cross_int[p3][ncross[p3]-1];
	    
	    ncross[p3]-=1;
	    ncrosscom[bases[p3].mycomplex]-=1;
	    j-=1;
	  }
	}
	
      }
    }
      
  }  
  
  ncross[p2]=-1;//this protein reacted
  ncross[p1]=-1;//this protein reacted
  ncrosscom[bases[p1].mycomplex]=-1;
  ncrosscom[bases[p2].mycomplex]=-1;
  

}
