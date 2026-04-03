#include <iostream>
#include <cmath>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include <gsl/gsl_sf_bessel.h>

//void transversePush_3D(Domain *D,int iteration);
void push_theta_gamma_1D(Domain &D,int iteration);
void push_theta_gamma_3D(Domain &D,int iteration);
//void drift_theta_gamma_3D(Domain *D,int iteration);
//void drift_theta_gamma_1D(Domain *D,int iteration);



void push_theta_gamma(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
    push_theta_gamma_1D(*D,iteration);
    break;

  case 2 :
//    particlePush2D(&D,iteration);
    break;
  case 3 :
    push_theta_gamma_3D(*D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}


void transversePush(Domain *D,int iteration)
{
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   int startI=1;  
   int endI=D->subSliceN+1;
   
   double dz=D->dz*0.5;
   double ku=D->ku;
   double K0=D->K0;
   double ue=D->ue;
   double aa=D->K0*D->K0*0.25*ue*D->ku;
   double bb=D->K0*D->K0*0.25*(1.0+ue*ue)*D->ku*D->ku;
   double cc=eCharge*D->g/eMass/velocityC;
   double xCoef=0.0;
   double yCoef=0.0;
   int K0_alpha=D->K0_alpha;

   if (D->undType==UndMode::BiPolar) {
      if (D->K0_alpha==1) {
         xCoef=1;
         yCoef=0;
      } else if (D->K0_alpha==-1) {
         xCoef=0;
         yCoef=1;
      } 
   } else if (D->undType==UndMode::QuadPolar) {
      xCoef=1;
      yCoef=1;
   } else {
      if(myrank==0) 
         std::cerr << "Error: No undulator_mode = " << K0_alpha << std::endl;
      MPI_Abort(MPI_COMM_WORLD,1);
   }
   
   double coefList[5]={0,0,0.5,0.5,1};
   double k_xL[5]={0,0,0,0,0};
   double k_yL[5]={0,0,0,0,0};
   double k_pxL[5]={0,0,0,0,0};
   double k_pyL[5]={0,0,0,0,0};

   for(int s=0; s<D->nSpecies; ++s)
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI) 
      {
         auto& p = D->particle[sliceI].head[s]->pt;         
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double gam=p->gamma[n];
            double x0=p->x[n];  
            double y0=p->y[n];
            double px0=p->px[n]; 
            double py0=p->py[n];
            double dd=ku*K0*K0/(4.0*gam*gam)*ue;
         
            k_xL[0]=0.0;
            k_yL[0]=0.0;
            k_pxL[0]=0.0;
            k_pyL[0]=0.0;
            
            for(int m=1; m<5; ++m) {
               double x = x0 + dz*k_xL[m-1]*coefList[m];
               double y = y0 + dz*k_yL[m-1]*coefList[m];
               double px = px0 + dz*k_pxL[m-1]*coefList[m];
               double py = py0 + dz*k_pyL[m-1]*coefList[m];

               k_xL[m]=px/gam-aa*(K0_alpha*x+y)/(gam*gam);
               k_yL[m]=py/gam+aa*(K0_alpha*x-y)/(gam*gam);
               k_pxL[m]=K0_alpha*(px-py)*dd-(bb/gam*xCoef+cc)*x;
               k_pyL[m]=K0_alpha*(px+py)*dd-(bb/gam*yCoef-cc)*y;
            }

            p->x[n]+=dz/6.0*(k_xL[1]+2*k_xL[2]+2*k_xL[3]+k_xL[4]);
            p->y[n]+=dz/6.0*(k_yL[1]+2*k_yL[2]+2*k_yL[3]+k_yL[4]);
            p->px[n]+=dz/6.0*(k_pxL[1]+2*k_pxL[2]+2*k_pxL[3]+k_pxL[4]);
            p->py[n]+=dz/6.0*(k_pyL[1]+2*k_pyL[2]+2*k_pyL[3]+k_pyL[4]);
            

         }
      }		//End of for(sliceI)
   }

}

void push_theta_gamma_3D(Domain &D,int iteration)
{
   ptclList *p;
   
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   
   int startI=1;       
   int endI=D.subSliceN+1;
   int numHarmony=D.numHarmony;
   int minI=D.minI;
   int maxI=D.maxI;
   int nx = D.nx;
   int ny = D.ny;
   int N=nx*ny;
   int L=D.SCLmode;

   double dz=D.dz;    
   double dx=D.dx; 
   double dy=D.dy; 
   double dr=D.dr;
   double ku=D.ku;    
   double ks=D.ks;
   double K0=D.K0;
   double dBessel = D.dBessel;
   double minX=D.minX;  
   double minY=D.minY;
   double dPhi=2*M_PI*D.numSlice;
   double e_mc2 = eCharge/eMass/velocityC/velocityC;	 

   double ue=D.ue;
   std::complex<double> U=(1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);
 
   double K0_alpha=D.K0_alpha;
   std::complex<double> etaX,etaY;
   if(K0_alpha==1) {
      etaX=std::cos(Phi*0.5);
      etaY=I*std::sin(Phi*0.5);
   } else if(K0_alpha==-1) {
      etaX=I*std::sin(Phi*0.5);
      etaY=std::cos(Phi*0.5);
   } else {
      if(myrank==0)
         std::cerr << "Error: Invalid K0_alpha = " << K0_alpha << std::endl;
      MPI_Abort(MPI_COMM_WORLD,1);
   }

   double coefList[5]={0,0,0.5,0.5,1},k_th[5]={0,0,0,0,0},k_gam[5]={0,0,0,0,0};	 
   double w[2]={0.0,0.0}, wx[2]={0.0,0.0}, wy[2]={0.0,0.0};
   double sumEzPart=0.0;
   //double xi = ks/ku*K0*K0*(1.0-ue*ue)/(8.0*gam*gam);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);

   cplx fx=0.0+I*0.0;
   cplx fy=0.0+I*0.0;
   std::vector<cplx> Ux(numHarmony),Uy(numHarmony),Em(L);

   for(int s=0; s<D.nSpecies; ++s)
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         p = D.particle[sliceI].head[s]->pt;         
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double x0=p->x[n];   
            double y0=p->y[n];
            double r2= x0*x0 + y0*y0;
            double px=p->px[n]; 
            double py=p->py[n];
            double pr2= px*px + py*py;
    	      double K2=1.0 + pr2 + (1.0+ue)*K0*K0*0.5*(1.0+ku*ku*0.5*r2);
            double th0=p->theta[n]; 
            double gam0=p->gamma[n];

            double xi2 = std::sqrt(2.0) * K0 / (1.0+K0*K0*0.5*(1.0+ue*ue))
                       * std::sqrt((1+K0_alpha)*(px*px+ue*ue*py*py)+(1-K0_alpha)*(px*px*ue*ue+py*py));
            double B=std::atan(
               2.0*K0_alpha*ue*(px-py)
               / ((1.0+K0_alpha)*(px*px+ue*ue*py*py)+(1-K0_alpha)*(px*px*ue*ue+py*py))
            );
            cplx expP=std::exp(I*(B-Phi*0.5));
            cplx expM=std::conj(expP);

            int idxI=(int)((x0-minX)/dx);
	         int idxJ=(int)((y0-minY)/dy);
            wx[1]=(x0-minX)/dx-idxI; wx[0]=1.0-wx[1];
            wy[1]=(y0-minY)/dy-idxJ; wy[0]=1.0-wy[1];
	         if(idxI>=0 && idxI<nx-1 && idxJ>=0 && idxJ<ny-1)  {
               /*
               double r=std::sqrt(r2);
               int idxR = r/dr;
               wr[1]=(r/dr-idxR); wr[0]=1.0-wr[1];
               if(idxR+1<D->nr) {
                  for(int ll=0; ll<L; ++ll) {
  		               Em[ll]=0.0+I*0.0;
                     for(f=0; f<F; f++) 
                        for(jj=0; jj<2; jj++) 
		                     Em[ll]+=D->Ez[sliceI][indexJ+jj][ll][f]*wr[jj];
                  }
		         }
               */       	   
               for(int h=0; h<numHarmony; ++h)  {				
                  Ux[h]=0.0+I*0.0;
                  Uy[h]=0.0+I*0.0;
                  for(int ii=0; ii<2; ++ii) 
                     for(int jj=0; jj<2; ++jj)  {
      			         Ux[h]+=D.Ux[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]*wx[ii]*wy[jj];
      			         Uy[h]+=D.Uy[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]*wx[ii]*wy[jj];
                     }
               }
               double sumU2=0.0;
               for(int h=0; h<numHarmony; ++h)  {
                  int H = D.harmony[h];
                  double absUx2 = std::norm(Ux[h]);
                  double absUy2 = std::norm(Uy[h]);
                  sumU2 += (absUx2 + absUy2)/(2.0*H*H);
               }
                 
               k_th[0]=0; 
               k_gam[0]=0; 
               for(size_t m=1; m<5; ++m) {
                  double th=th0 + dz*k_th[m-1]*coefList[m]; 
                  double gam=gam0 + dz*k_gam[m-1]*coefList[m];
                     
                  double sumTh=0.0;
                  double sumG=0.0;

                  sumEzPart = 0.0;
                  //for(int ll=0; ll<L; ++ll)  {
                  //   double tmp=std::real(Em[ll]*std::exp(I*(ll+1.0)*th));
                  //   sumEzPart += 2.0*tmp;
                  //}

                  for(int h=0; h<numHarmony; ++h)  {
                     int H = D.harmony[h];
                     int idx=H*xi/dBessel;
                     w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                     if(H%2==1)  {  //odd harmony
                        double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
                        int order=(H-1)/2;
                        double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        order=(H+1)/2;
                        double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        fx=sign*(K0_alpha*J1-J2);
                        fy=sign*(K0_alpha*J1+J2);
                     } else {    //even harmony
                        double sign = (H/2 % 2 ==0) ? 1.0 : -1.0;
                        int idx=(H*xi)/dBessel;
                        w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                        int order=(H-2)/2;
                        double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        order=H/2;
                        double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        order=(H+2)/2;
                        double J3=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                        fx=sign*xi2*H*0.5*(expP*(K0_alpha*J1-J2)+expM*(K0_alpha*J2-J3));
                        fy=sign*xi2*H*0.5*(expP*(K0_alpha*J1+J2)+expM*(K0_alpha*J2+J3));
                     }
                     cplx tmpComp=std::exp(I*(1.0*H)*(th+Phi*0.5))*(etaX*fx*Ux[h]-etaY*fy*Uy[h]);
                     sumTh += std::real(I*tmpComp)/(1.0*H);
                     sumG += std::real(tmpComp)/(1.0*gam);
                  }  //End of harmonics

                  k_th[m]=ku-ks/(2.0*gam*gam)*(K2+sumU2+K0*std::sqrt(ue*ue+1.0)*sumTh);
                  k_gam[m]=ks*K0*std::sqrt(ue*ue+1.0)*sumG + e_mc2*sumEzPart;
               }   //End of Runge-Kutta
               
               double tmpTh=dz/6.0*(k_th[1] + 2*k_th[2] + 2*k_th[3] + k_th[4]);
               if(tmpTh>dPhi || tmpTh<-dPhi) {
                  printf("myrank=%d,iteration=%d,dTheta=%g,sumEzPart=%g,i=%d,j=%d,k_th[1]=%g,k_th[2]=%g,k_th[3]=%g,k_th[4]=%g,tmpTh=%g\n",myrank,iteration,dPhi,sumEzPart,idxI,idxJ,k_th[1],k_th[2],k_th[3],k_th[4],tmpTh);
                     exit(0);
               } else;
               p->theta[n]+=tmpTh; 
               double tmpGam=dz/6.0*(k_gam[1]+2*k_gam[2]+2*k_gam[3]+k_gam[4]); //-dz*wakeE; 
               p->gamma[n]+=tmpGam;
            }  // End of if(idxI,idxJ)
         }   // End of for(cntN)

      }     //Enf of for(sliceI)     

   }        //Enf of for(s)

}




void push_theta_gamma_1D(Domain &D,int iteration)
{
   ptclList *p;

   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI=1;       
   int endI=D.subSliceN+1;
   int minI=D.minI;   
   int maxI=D.maxI;
   double ue=D.ue;
   double ks=D.ks;
   double ku=D.ku;
   double K0=D.K0;
   double dz=D.dz;
   int L = D.SCLmode;
   double dPhi=2.0*M_PI*D.numSlice;

   std::complex<double> U=(1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);
   int numHarmony = D.numHarmony;
   double e_mc2 = eCharge/eMass/velocityC/velocityC;	 

   double K0_alpha = D.K0_alpha;
   std::complex<double> etaX,etaY;
   if(K0_alpha==1) {
      etaX=std::cos(Phi*0.5);
      etaY=I*std::sin(Phi*0.5);
   } else if(K0_alpha==-1) {
      etaX=I*std::sin(Phi*0.5);
      etaY=std::cos(Phi*0.5);
   } else {
      if(myrank==0) 
         std::cerr << "Error: Invalid K0_alpha = " << K0_alpha << std::endl;
      MPI_Abort(MPI_COMM_WORLD,1);
   }

   
   double K2 = 1.0 + K0*K0*0.5*(1.0+ue*ue);
   double coefList[5]={0,0,0.5,0.5,1},k_th[5]={0,0,0,0,0},k_gam[5]={0,0,0,0,0},w[2]={0,0};
   double dBessel = D.dBessel;
   double sumEzPart=0.0;
   //double xi = ks/ku*K0*K0*(1.0-ue*ue)/(8.0*gam*gam);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);
   cplx fx=0.0+I*0.0;
   cplx fy=0.0+I*0.0;
   std::vector<cplx> Ux(numHarmony),Uy(numHarmony),Em(L);

   for(int s=0; s<D.nSpecies; ++s)
   {
      for(int sliceI=startI; sliceI<endI; ++sliceI) 
      {
         for(int ll=0; ll<L; ++ll) 
	         Em[ll]=D.Ez[ll][sliceI];
         
         for(int h=0; h<numHarmony; ++h)  {				
            Ux[h]=D.Ux[h][sliceI];
            Uy[h]=D.Uy[h][sliceI];
         }
               
         double sumU2=0.0;
         for(int h=0; h<numHarmony; ++h)  {
            int H = D.harmony[h];
            double absUx2 = std::norm(Ux[h]);
            double absUy2 = std::norm(Uy[h]);
            sumU2 += (absUx2 + absUy2)/(2.0*H*H);
         }
               
         p = D.particle[sliceI].head[s]->pt;         
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double th0=p->theta[n]; 
            double gam0=p->gamma[n];
            
            //Start of Runge-Kutta
            k_th[0]=0.0; 
            k_gam[0]=0.0; 
            for(size_t m=1; m<5; ++m) {
               double th=th0 + dz*k_th[m-1]*coefList[m]; 
               double gam=gam0 + dz*k_gam[m-1]*coefList[m];
               double sumTh=0.0;
               double sumG=0.0;
         
               sumEzPart = 0.0;
               for(int ll=0; ll<L; ++ll)  {
                  double tmp=std::real(Em[ll]*std::exp(I*(ll+1.0)*th));
                  sumEzPart += 2.0*tmp;
               }

               for(int h=0; h<numHarmony; ++h)  {
                  int H = D.harmony[h];
                  double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
                  int idx=(int)(H*xi/dBessel);
                  w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                  if(H%2==1)  {  //odd harmony
                     int order=(H-1)*0.5;
                     double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=(H+1)*0.5;
                     double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     fx=sign*(K0_alpha*J1-J2);
                     fy=sign*(K0_alpha*J1+J2);
                  } else {    //even harmony
                     fx=0.0+I*0.0;
                     fy=0.0+I*0.0;
                  }

                  cplx tmpComp=std::exp(I*(1.0*H)*(th+Phi*0.5))*(etaX*fx*Ux[h]-etaY*fy*Uy[h]);
                  sumTh += std::real(I*tmpComp)/(1.0*H);
                  sumG += std::real(tmpComp)/(1.0*gam);
               }  //End of harmonics

               k_th[m]=ku-ks/(2.0*gam*gam)*(K2+sumU2+K0*sqrt(ue*ue+1.0)*sumTh);
               k_gam[m]=ks*K0*std::sqrt(ue*ue+1.0)*sumG + e_mc2*sumEzPart;
            }   //End of Runge-Kutta 
   
            double tmpTh=dz/6.0 * (k_th[1]+2.0*k_th[2]+2.0*k_th[3]+k_th[4]);
            if(tmpTh>dPhi || tmpTh<=-dPhi) {
               printf("myrank=%d,iteration=%d,sliceI=%d,th0=%g,dTheta=%g,sumEzPart=%g,k_th[1]=%g,k_th[2]=%g,k_th[3]=%g,k_th[4]=%g\n"
                   ,myrank,iteration,sliceI,th0,dPhi,sumEzPart,k_th[1],k_th[2],k_th[3],k_th[4]);
               MPI_Abort(MPI_COMM_WORLD,1);
            }
            p->theta[n]+=tmpTh; 
            double tmpGam=dz/6.0*(k_gam[1]+2*k_gam[2]+2*k_gam[3]+k_gam[4]); //-dz*wakeE; 
            p->gamma[n]+=tmpGam;
         }
         
      }  //End of for(sliceI)

   }   //End of for(s)


}





/*

void dS_dz(double *S,double *dS,double *coef)
{
   double X=S[0], Y=S[1], px=S[2], py=S[3];
   double aa=coef[1], bb=coef[2], cc=coef[3];
   double invG=1.0/coef[0];
   dS[0] = px*invG - aa*(X+Y);
   dS[1] = py*invG + aa*(X-Y);
   dS[2] = (px-py)*aa - (bb+cc)*X;
   dS[3] = (px+py)*aa - (bb-cc)*Y;
}

void rk4_S(double *S, double dz,double *coef) {
   double k1[4], k2[4], k3[4], k4[4];
   double temp[4];
   int i;

   dS_dz(S, k1, coef);

   for (i = 0; i < 4; i++) temp[i] = S[i] + dz*0.5 * k1[i];
   dS_dz(temp, k2, coef);

   for (i = 0; i < 4; i++) temp[i] = S[i] + dz*0.5 * k2[i];
   dS_dz(temp, k3, coef);

   for (i = 0; i < 4; i++) temp[i] = S[i] + dz * k3[i];
   dS_dz(temp, k4, coef);

   for (i = 0; i < 4; i++) 
      S[i] += (dz/6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}



void transversePush(Domain *D,int iteration)
{
  switch(D->dimension)  {
  case 1 :
	  ;
//    particlePush1D(D);
    break;

  case 2 :
//    particlePush2D(&D,iteration);
    break;
  case 3 :
    transversePush_3D(D,iteration);
    break;
    
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}



void drift_theta_gamma(Domain *D,int iteration)
{
   int startI,endI,minI,maxI,sliceI,s,numInBeamlet,nx,ny,n,i,j,h,ii,jj,numHarmony;
   LoadList *LL;
   double dz,dx,dy,ks,ku,x0,y0,px,py,pr2,invGam,invGam0,wakeE,tmp,invBeta0;
   double minX,minY,wx[2],wy[2],sumU2;
   double complex Ux[D->numHarmony],Uy[D->numHarmony];
   ptclList *p;

   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   dz=D->dz;    ks=D->ks;   ku=D->ku;
   dx=D->dx; dy=D->dy;
   nx=D->nx; ny=D->ny;

   startI=1;       endI=D->subSliceN+1;
   numHarmony=D->numHarmony;
   minI=D->minI;   maxI=D->maxI;
   minX=D->minX;  minY=D->minY;

   invGam0=1.0/D->gamma0;
   invBeta0=1.0/sqrt(1-invGam0*invGam0);

   if(myrank==0) printf("iteration=%d, driftFlag=%d, drift ON, K0=%g\n",iteration,D->driftFlag,D->K0);

   LL=D->loadList;
   s=0;
   while(LL->next) {
      numInBeamlet=LL->numInBeamlet;
      for(sliceI=startI; sliceI<endI; sliceI++)
      {
         if(D->wakeONOFF==ON) wakeE=D->wakeE[sliceI-startI+minI]/mc2*1e-6;
         else                 wakeE=0.0;

         p=D->particle[sliceI].head[s]->pt;
         while(p) {
            for(n=0; n<numInBeamlet; n++) {
               px=p->px[n];  py=p->py[n];  invGam=1.0/p->gamma[n];
               //tmp=ks*0.5*dz*invBeta0*(invGam0*invGam0-invGam*invGam*(1+px*px+py*py));
               x0=p->x[n];   y0=p->y[n];
               px=p->px[n]; py=p->py[n];
               pr2= px*px + py*py;

               i=(int)((x0-minX)/dx);
               j=(int)((y0-minY)/dy);
               wx[1]=(x0-minX)/dx-i; wx[0]=1.0-wx[1];
               wy[1]=(y0-minY)/dy-j; wy[0]=1.0-wy[1];
               sumU2=0.0;
	            if(i>=0 && i<nx-1 && j>=0 && j<ny-1)  {
                  for(h=0; h<numHarmony; h++)  {				
                     Ux[h]=0.0+I*0.0;
                     Uy[h]=0.0+I*0.0;
                     for(ii=0; ii<2; ii++) 
                        for(jj=0; jj<2; jj++)  {
      			            Ux[h]+=D->Ux[h][sliceI][(j+jj)*nx+(i+ii)]*wx[ii]*wy[jj];
      			            Uy[h]+=D->Uy[h][sliceI][(j+jj)*nx+(i+ii)]*wx[ii]*wy[jj];
                        }
                  }
                  for(h=0; h<numHarmony; h++)
                     sumU2+=cabs(Ux[h])*cabs(Ux[h])+cabs(Uy[h])*cabs(Uy[h]);
               }
               tmp=ku-ks*0.5*invGam*invGam*(1+pr2+sumU2);
               p->theta[n]+=tmp*dz;
               p->gamma[n]-=wakeE*dz;
            }
            p=p->next;
         }   //Enf of while(p)			
      }     //End of for(sliceI)
		
      LL=LL->next;
      s++;
   }
}





/*
void push_theta_gamma_1D(Domain *D,int iteration)
{
   int numHarmony,order,startI,endI,minI,maxI;
   int n,i,s,h,H,ll,L,idx,intThe,bn,m,numInBeamlet;	
   LoadList *LL;
   double complex Ux[D->numHarmony],Uy[D->numHarmony],Em[D->SCLmode],compVal;
   double dz,ku,ks,K0,K,xi,e_mc2,dBessel,w[2];
   double z,gamma,theta,invGam,invBeta,tmp,dPhi,sumGam,sumTh,sumEzPart,prevThe;
	double coef,lkCoef,JJ,J1,J2,wakeE,kList[5],lList[5];
   ptclList *p;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   L = D->SCLmode;
   startI=1;       endI=D->subSliceN+1;
   minI=D->minI;   maxI=D->maxI;
	dz = D->dz;     K0=D->K0;
   ku=D->ku;       ks=D->ks;
   numHarmony = D->numHarmony;
   dBessel = D->dBessel;
   dPhi=2*M_PI*D->numSlice;
   e_mc2 = eCharge/eMass/velocityC/velocityC;	 
	bn=D->bn;

   LL=D->loadList;
   s=0;
   while(LL->next) {
      numInBeamlet=LL->numInBeamlet;

      for(i=startI; i<endI; i++)
      {
         if(D->wakeONOFF==ON) wakeE=D->wakeE[i-startI+minI]/mc2*1e-6;
         else                wakeE=0.0;      

         for(ll=0; ll<L; ll++) Em[ll]=D->Ez[i][0][ll][0];
         for(h=0; h<numHarmony; h++) {
            Ux[h]=D->Ux[h][i][0];
            Uy[h]=D->Uy[h][i][0];
         }

         K=K0;
         p=D->particle[i].head[s]->pt;
         while(p) {
            for(n=0; n<numInBeamlet; n++) {
	            kList[0]=0.0;
		         lList[0]=0.0;
               for(m=1; m<5; m++) {
                  lkCoef=((int)(m*0.5))*0.5;
	               theta=p->theta[n] + lkCoef*kList[m-1]*dz;
	               gamma=p->gamma[n] + lkCoef*lList[m-1]*dz; invGam=1.0/gamma;
	               invBeta = 1.0-(1.0 + K*K*0.5)*0.5*invGam*invGam;
	               invBeta = 1.0/invBeta;
	               sumTh=sumGam=0.0;
	               //xi=ks/ku*0.25*K*K*invGam*invGam;
	               xi=K*K*0.25/(1+K*K*0.5);          
			         for(h=0; h<numHarmony; h++)  {
                     H = D->harmony[h];
                     if(H%2==1) {
                        coef=pow(-1.0,(H-1)*0.5);
                        idx=(int)(H*xi/dBessel);
                        if(idx>bn-1) { printf("idx=%d\n",idx); idx=bn-2; }
                        w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                        order=(H-1)*0.5;
                        J1=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                        order=(H+1)*0.5;
                        J2=D->BesselJ[idx][order]*w[0]+D->BesselJ[idx+1][order]*w[1];
                        JJ=coef*(J1-J2);
	                  } else JJ=0.0;
                     compVal=Ux[h]*cexp(I*H*theta);
                     sumTh -=2*JJ*K/sqrt(2.0)*cimag(compVal);
		     sumGam-=2*JJ*K/sqrt(2.0)*creal(compVal);
                  }
                  sumEzPart = 0.0;
                  for(ll=0; ll<L; ll++)  {
                     tmp=creal(Em[ll]*cexp(I*(ll+1)*theta));
                     sumEzPart += 2.0*tmp;
                  }
                  kList[m]=dz*(ku-ks*(1.0+K*K*0.5+sumTh)*0.5*invGam*invGam)*invBeta;
						lList[m]=dz*(ks*sumGam/2.0*invGam*invBeta + e_mc2*sumEzPart);
			      }  //End of for(m)

               tmp=dz/6.0*(kList[1]+2*kList[2]+2*kList[3]+kList[4]);
               if(tmp>dPhi || tmp<-dPhi) {
                  printf("iteration=%d, dTheta=%g, sumEzPart=%g, Ux[0]=%g\n",iteration,tmp,sumEzPart,cabs(D->Ux[0][i][0]));  //lala
                  exit(0);
               } else ;
               p->theta[n]+=tmp;
               p->gamma[n]+=dz/6.0*(lList[1]+2*lList[2]+2*lList[3]+lList[4]) - wakeE*dz;
		      }
	         p=p->next;
         }    //End of while9p)
      }		//End of for(i)
		LL=LL->next;
		s++;
   }
}
*/
//lala
/*	
void push_theta_gamma_1D(Domain *D,int iteration)
{
    int i,s,h,harmony,order,startI,endI,minI,maxI;
    LoadList *LL;
    double complex U[D->harmony],Ez[D->harmony],compValU,compValEz,expITh,expIhTh;
    double dz,ku,ks,K0,K,xi,ar[D->harmony],psi[D->harmony];
    double z,gamma,theta,invGam,invBeta,tmp,w,dPhi;
    double coef,JJ[D->harmony],J1,J2,sumTh,sumGam,wakeE;
    double k1,k2,k3,k4,l1,l2,l3,l4,totalEz;
    ptclList *p;

    startI=1;       endI=D->subSliceN+1;
    minI=D->minI;   maxI=D->maxI;
    dz=D->dz;
    ku=D->ku;
    ks=D->ks;
    K0=D->K0;
    harmony=D->harmony;
    dPhi=2*M_PI*D->numSlice;

    for(i=startI; i<endI; i++)
    {
      for(h=0; h<harmony; h++) { 
        ar[h]=cabs(D->U[h][i][0]);
        psi[h]=carg(I*D->U[h][i][0]);
        Ez[h]=D->Ez[h][i][0];
      }
      if(D->wakeONOFF=ON) wakeE=D->wakeE[i-startI+minI]/mc2*1e-6;
      else                wakeE=0.0;      

      for(s=0; s<D->nSpecies; s++)  {
        p=D->particle[i].head[s]->pt;
        while(p) {
          K=K0;
	  xi=K*K*0.5/(1+K*K);

          for(h=1; h<=harmony; h++)  {
            if(h%2==1)  {  //odd harmony
              coef=pow(-1.0,(h-1)/2);
              order=(h-1)/2;
              J1=gsl_sf_bessel_Jn(order,h*xi);
              order=(h+1)/2;
              J2=gsl_sf_bessel_Jn(order,h*xi);
              JJ[h-1]=coef*(J1-J2);
            } else {
              JJ[h-1]=0.0;
            }
          }

	  theta=p->theta;
	  gamma=p->gamma; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*2*ar[0]*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;	  
          k1=ku+ks*(1.0-invBeta);
	  l1=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

	  theta=p->theta+0.5*dz*k1;
	  gamma=p->gamma+0.5*dz*l1; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*ar[0]*2*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;
          k2=ku+ks*(1.0-invBeta);
	  l2=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

	  theta=p->theta+0.5*dz*k2;
	  gamma=p->gamma+0.5*dz*l2; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*ar[0]*2*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;
          k3=ku+ks*(1.0-invBeta);
	  l3=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

	  theta=p->theta+dz*k3;
	  gamma=p->gamma+dz*l3; invGam=1.0/gamma;
//	  invBeta=1.0+(1+K*K)*0.5*invGam*invGam;
	  invBeta=1.0+(1+K*K+K*JJ[0]*ar[0]*2*cos(theta+psi[0]))*0.5*invGam*invGam;
          compValEz=0.0+I*0.0;
          for(h=1; h<=harmony; h++)  {
	    expIhTh=cexp(I*h*theta);
            compValEz+=Ez[h-1]*expIhTh;
	  }
	  totalEz=eCharge/eMass/velocityC/velocityC*2*creal(compValEz)-wakeE;
          k4=ku+ks*(1.0-invBeta);
	  l4=-ks*JJ[0]*K*ar[0]*invBeta*invGam*sin(theta+psi[0])+totalEz;

          tmp=dz*1.0/6.0*(k1+2*k2+2*k3+k4);
	  if(tmp>dPhi || tmp<-1.0*dPhi) {
		  printf("tmp=%g,dPhi=%g\n",tmp,dPhi);
	  } else ;
          p->theta+=tmp;
          p->gamma+=dz*1.0/6.0*(l1+2*l2+2*l3+l4);

//          // energy loss by wake field
//          delGam=D->wakeE[i-startI+minI]/mc2;	  
//          p->gamma-=delGam*dz;

          tmp=p->theta/dPhi;
	  w=tmp-(int)tmp;
	  p->theta=w*dPhi;
          if(p->theta>dPhi)   p->theta=p->theta-dPhi;
          else if(p->theta<0) p->theta=dPhi+p->theta;
	  

	  p=p->next;
        }
      }		//End of for(s)
    }		//End of for(i)
}
*/


/*
double Runge_Kutta_gamma(double complex *Ux,double theta,double harmony,double ks,double K,double xi,double gamma,double invBeta0,double dz)
{
  double sum,sinX,cosX,coef,J1,J2,JJ,k1,k2,k3,k4;
  double complex compVal;
  int h,order;

  sum=0.0;
  sinX=sin(theta); cosX=cos(theta);
  for(h=1; h<=harmony; h++)  {
    if(h%2==1)  {  //odd harmony
      coef=pow(-1.0,h-1);
      order=(h-1)/2;
      J1=gsl_sf_bessel_Jn(order,h*xi);
      order=(h+1)/2;
      J2=gsl_sf_bessel_Jn(order,h*xi);
      JJ=coef*(J1-J2);
      compVal=Ux[h-1]*cexp(-I*theta);
      sum=JJ*2*creal(compVal)*h;
    } else {
      JJ=0.0;
      compVal=Ux[h-1]*cexp(-I*theta);
      sum=JJ*2*creal(compVal)*h;
    }
  }

  k1=-1.0*K*sqrt(2.0)*ks*invBeta0*sum/(gamma);
  
  k2=-1.0*K*sqrt(2.0)*ks*invBeta0*sum/(gamma+0.5*k1);
  
  k3=-1.0*K*sqrt(2.0)*ks*invBeta0*sum/(gamma+0.5*k2);
  
  k4=-1.0*K*sqrt(2.0)*ks*invBeta0*sum/(gamma+k3);

  return (k1/6.0+k2/3.0+k3/3.0+k4/6.0)*dz;
}

void phaseShift(Domain *D,int iteration)
{
   int n,s,i,sliceI,startI,endI,numInBeamlet;
   double shiftValue,theta;
   LoadList *LL;
   ptclList *p;
   PhaseShifter *PS;
   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks); 

   startI=1;       endI=D->subSliceN+1;

   PS=D->psList;
   while(PS->next) {
      for(i=0; i<PS->num; i++) {
         if(iteration==PS->step[i]) {
            shiftValue=PS->phase;
            if(myrank==0) printf("phase shift with %g is done at step%d.\n",shiftValue,iteration);  else ;

            LL=D->loadList;
            s=0;
            while(LL->next) {
               numInBeamlet=LL->numInBeamlet;
        
               for(sliceI=startI; sliceI<endI; sliceI++)  {
                  p=D->particle[sliceI].head[s]->pt;
                  while(p) {
                     for(n=0; n<numInBeamlet; n++) p->theta[n]-=shiftValue;
                     p=p->next;
                  }
               }		//End of for(sliceI)

               LL=LL->next;
               s++;
            }
         } else ;
      }
      PS=PS->next;
   }

}
*/

