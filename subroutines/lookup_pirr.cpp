#include "Vrnumer.h"
#include "numeric_GF.h"
#include "GF_calls.h"

double lookup_pirr(double R1, double r0, double deltat, double rc, double Dtot, double bindrad, double **FitMat, Vrnumer &limit, double passoc )
{
  double pirr, pirr1;
  double alpha;
  /*Establish which limits to use based on r0 and r1. 
    need to calc closest grid points for r0, and interpolate.
    
  */
  int r0_downind, r0_upind;
  double r0_down, r0_up;
  double Dfit, kfit, Sfit, rcfit;
  double slope;
  double delr0=limit.delr0;
  double pirrdown, pirrup;
  double kdiff;
  r0_downind=int((r0-bindrad)/delr0);
  r0_upind=r0_downind+1;
  r0_down=r0_downind*delr0+bindrad;
  r0_up=r0_upind*delr0+bindrad;
  
  if(r0<limit.r0_2hybridstart){
    /*perform double fit to pirrev*/
    if(R1<FitMat[r0_downind][6]){
      /*use first half D, k, and S*/
      
      Dfit=FitMat[r0_downind][0];
      kfit=FitMat[r0_downind][1];
      Sfit=FitMat[r0_downind][2];
      
      kdiff=4*M_PI*Dfit*bindrad;
      alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
      pirr1= pirrev_valueF(R1, r0_down, deltat,  Dfit,  bindrad,  alpha);
      pirrdown=pirr1*Sfit;
      
      Dfit=FitMat[r0_upind][0];
      kfit=FitMat[r0_upind][1];
      Sfit=FitMat[r0_upind][2];
            
      kdiff=4*M_PI*Dfit*bindrad;
      alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
      pirr1= pirrev_valueF(R1, r0_up,deltat,  Dfit,  bindrad,  alpha);
      pirrup=pirr1*Sfit;

      slope=(pirrup-pirrdown)/delr0;
      pirr=pirrdown+slope*(r0-r0_down);
      //cout <<"In first regime, r0: "<<r0<<" 1st half R1: "<<R1<<" r0low, index: "<<r0_downind<<" r0_lowd: "<<r0_down<<" up: "<<r0_upind<<'\t';
      //cout <<"pirrup: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;
    }else{
      /*use second half D, k, S in 2pirr regime*/
      
  
      Dfit=FitMat[r0_downind][3];
      kfit=FitMat[r0_downind][4];
      Sfit=FitMat[r0_downind][5];
      // cout <<"In first regime, r0: "<<r0<<" second half R1: "<<R1<<" r0low, index: "<<r0_downind<<" r0_lowd: "<<r0_down<<" up: "<<r0_upind<<'\t';
      kdiff=4*M_PI*Dfit*bindrad;
      alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
      pirr1= pirrev_valueF(R1, r0_down, deltat,  Dfit,  bindrad,  alpha);
      pirrdown=pirr1*Sfit;
      if(r0+delr0>=limit.r0_2hybridstart){
	/*upper bin is in hybrid regime, use free fit*/
	Dfit=FitMat[r0_upind][3];
	rcfit=FitMat[r0_upind][4];
	Sfit=FitMat[r0_upind][5];
	pirr1=vr_pfree_value_norm(R1, r0_up, deltat, Dfit, bindrad, rcfit);
	pirrup=pirr1*Sfit;
	slope=(pirrup-pirrdown)/delr0;
	pirr=pirrdown+slope*(r0-r0_down);
	//cout <<"pirrup, in free range!: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;
      }else{
	//use pirr fit
	Dfit=FitMat[r0_upind][3];
	kfit=FitMat[r0_upind][4];
	Sfit=FitMat[r0_upind][5];
	
	
	kdiff=4*M_PI*Dfit*bindrad;
	alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
	pirr1= pirrev_valueF(R1, r0_up, deltat,  Dfit,  bindrad,  alpha);
	pirrup=pirr1*Sfit;
	slope=(pirrup-pirrdown)/delr0;
	pirr=pirrdown+slope*(r0-r0_down);
	//cout <<"pirrup: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;
      }
      
    }
  }else if(r0<limit.r0_2freefitstart){
    /*Second regime, first half fits to pirr, second half to pfree*/
    if(R1<FitMat[r0_downind][6]){
      
      Dfit=FitMat[r0_downind][0];
      kfit=FitMat[r0_downind][1];
      Sfit=FitMat[r0_downind][2];
      //cout <<"In second regime, r0: "<<r0<<" first half R1: "<<R1<<" r0low, index: "<<r0_downind<<" r0_lowd: "<<r0_down<<" up: "<<r0_upind<<'\t';
      kdiff=4*M_PI*Dfit*bindrad;
      alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
      pirr1= pirrev_valueF(R1, r0_down, deltat,  Dfit,  bindrad,  alpha);
      pirrdown=pirr1*Sfit;
      
      if(r0+delr0>limit.r0_2freefitstart){
	/*upper bin switches to pfree*/
	Dfit=FitMat[r0_upind][0];
	rcfit=FitMat[r0_upind][1];
	Sfit=FitMat[r0_upind][2];
	pirr1=vr_pfree_value_norm(R1, r0_up, deltat, Dfit, bindrad, rcfit);
	pirrup=pirr1*Sfit;
	slope=(pirrup-pirrdown)/delr0;
	pirr=pirrdown+slope*(r0-r0_down);
	//cout <<"pirrup in free zone!: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;
      }else{
	Dfit=FitMat[r0_upind][0];
	kfit=FitMat[r0_upind][1];
	Sfit=FitMat[r0_upind][2];
	kdiff=4*M_PI*Dfit*bindrad;
	alpha=(1+kfit/kdiff)*sqrt(Dfit)/bindrad;
	pirr1= pirrev_valueF(R1, r0_up,deltat,  Dfit,  bindrad,  alpha);
	pirrup=pirr1*Sfit;
	slope=(pirrup-pirrdown)/delr0;
	pirr=pirrdown+slope*(r0-r0_down);
	//	cout <<"pirrup: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;
      }
    }else{
      /*second half hybrid, Use pfree fit, with scale fact*/
      
      Dfit=FitMat[r0_downind][3];
      rcfit=FitMat[r0_downind][4];
      Sfit=FitMat[r0_downind][5];
      pirr1=vr_pfree_value_norm(R1, r0_down, deltat, Dfit, bindrad, rcfit);
      pirrdown=pirr1*Sfit;
      //cout <<"In second regime, 2nd half  pirrdown in free "<<pirrdown<<endl;
    
      Dfit=FitMat[r0_upind][3];
      rcfit=FitMat[r0_upind][4];
      Sfit=FitMat[r0_upind][5];
      //cout <<"pirr up in freefit regime, r0: "<<r0<<" R1: "<<R1<<" r0low, index: "<<r0_downind<<" r0_lowd: "<<r0_down<<" up: "<<r0_upind<<'\t';
      pirr1=vr_pfree_value_norm(R1, r0_up, deltat, Dfit, bindrad, rcfit);
      pirrup=pirr1*Sfit;
      slope=(pirrup-pirrdown)/delr0;
      pirr=pirrdown+slope*(r0-r0_down);
      //cout <<"pirrup: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;
    }
  }else if(r0<limit.r0_freeonlystart){
    /*fit parms of pirr are derived for pfree*/
    
    if(R1<FitMat[r0_downind][6]){
      /*bin below use free fit.*/
      Dfit=FitMat[r0_downind][0];
      rcfit=FitMat[r0_downind][1];
      Sfit=FitMat[r0_downind][2];
      pirr1=vr_pfree_value_norm(R1, r0_down, deltat, Dfit, bindrad, rcfit);
      pirrdown=pirr1*Sfit;
      
      /*if upper bin is greater than limit, we still have parameters for itt*/
      Dfit=FitMat[r0_upind][0];
      rcfit=FitMat[r0_upind][1];
      Sfit=FitMat[r0_upind][2];
      //cout <<"In 2freefit (3rd) regime, r0: "<<r0<<" first half R1: "<<R1<<" r0low, index: "<<r0_downind<<" r0_lowd: "<<r0_down<<" up: "<<r0_upind<<'\t';
      pirr1=vr_pfree_value_norm(R1, r0_up, deltat, Dfit, bindrad, rcfit);
      pirrup=pirr1*Sfit;
      slope=(pirrup-pirrdown)/delr0;
      pirr=pirrdown+slope*(r0-r0_down);
      //cout <<"pirrup: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;
    }else{
      /*Second half, use free fit up or down*/
      Dfit=FitMat[r0_downind][3];
      rcfit=FitMat[r0_downind][4];
      Sfit=FitMat[r0_downind][5];
      pirr1=vr_pfree_value_norm(R1, r0_down, deltat, Dfit, bindrad, rcfit);
      pirrdown=pirr1*Sfit;
      Dfit=FitMat[r0_upind][3];
      rcfit=FitMat[r0_upind][4];
      Sfit=FitMat[r0_upind][5];
      //cout <<"In 2freefit (3rd) regime, r0: "<<r0<<" second half R1: "<<R1<<" r0low, index: "<<r0_downind<<" r0_lowd: "<<r0_down<<" up: "<<r0_upind<<'\t';
      pirr1=vr_pfree_value_norm(R1, r0_up, deltat, Dfit, bindrad, rcfit);
      pirrup=pirr1*Sfit;
      slope=(pirrup-pirrdown)/delr0;
      pirr=pirrdown+slope*(r0-r0_down);
      //cout <<"pirrup: "<<pirrup<<" pirrdown "<<pirrdown<<" sampled: "<<pirr<<endl;

    }
  }else{
    /*if r0 gets large, we can't distinguish numerical noise, so set pirr=pfree_vr, Ratio=1*/
    pirr=vr_pfree_value_norm(R1, r0, deltat, Dtot, bindrad, rc);
    pirr*=(1.0-passoc);
    //cout <<"Setting pirr to pfree_vr: "<<pirr<<" r0: "<<r0<<" R1: "<<R1<< endl;
  }
  return pirr;

}
