#include <iostream>
#include <cmath>
#include "mesh.h"
#include "constants.h"
#include <gsl/gsl_sf_bessel.h>
#include <mpi.h>

//void solveField1D(Domain *D,int sliceI,int iteration);
void solve_Sc_1D(Domain &D,int iteration);
void solve_Field_U_1D(Domain &D,int iteration);
void solve_Field_U_3D(Domain &D,
                      std::vector<std::vector<cplx>>& Un,
                      std::vector<std::vector<cplx>>& Sc,
                      int iteration);
void solve_Sc_3D(Domain &D,int iteration);
//void MPI_Transfer1F_Z(double complex ***f1,int harmony,int N,int fromI,int toI);
//void MPI_Transfer1F_Zplus(double complex ***f1,int harmony,int N,int fromI,int toI);

/*
void shiftField(Domain D,int iteration)
{
   int h,numHarmony,i,j,startI,endI,N;

   N=D.nx*D.ny;
   numHarmony=D.numHarmony;
   startI=1;  endI=D.subSliceN+1;

   int myrank, nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   for(h=0; h<numHarmony; h++)  
      for(i=endI; i>=startI; i--) 
         for(j=0; j<N; j++) {
            D.Ux[h][i][j]=D.Ux[h][i-1][j];
            D.Uy[h][i][j]=D.Uy[h][i-1][j];
         }
		   
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Transfer1F_Zplus(D.Ux,D.numHarmony,N,endI,startI);
   MPI_Transfer1F_Zplus(D.Uy,D.numHarmony,N,endI,startI);
}
*/



void solveField(Domain *D,int iteration)
{
   int startI = 0;
   int endI = D->subSliceN+2;

   for(int h=0; h<D->numHarmony; ++h)
      for(int sliceI=0; sliceI<endI; ++sliceI) {
         D->ScUx[h][sliceI]=0.0+I*0.0;
         D->ScUy[h][sliceI]=0.0+I*0.0;
      }

   switch(D->dimension) {

   //1D field
   case 1:
      if(D->currentFlag==true) 
         solve_Sc_1D(*D,iteration);
      solve_Field_U_1D(*D,iteration);
      break;
   case 3:
      if(D->currentFlag==true) 
         solve_Sc_3D(*D,iteration);
      solve_Field_U_3D(*D,D->Ux,D->ScUx,iteration);
      solve_Field_U_3D(*D,D->Uy,D->ScUy,iteration);
      break;

   default:
      printf("In EzSolve.c, what dimension?\n");
   }
}

void solve_Field_U_3D(Domain &D,
                      std::vector<std::vector<cplx>>& Un,
                      std::vector<std::vector<cplx>>& Sc,
                      int iteration)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI=1;   
   int endI=D.subSliceN+1;
   int numHarmony=D.numHarmony;
   int nx=D.nx;  
   int ny=D.ny;
   int N=nx*ny;

   double dx=D.dx;  
   double dy=D.dy;
   double dz=D.dz;
   double ks=D.ks;
   double  currentFlag=D.currentFlag ? 1.0 : 0.0;
   int Lx=D.abcN; 
   int Ly=D.abcN;
   double sig0=D.abcSig;

   std::vector<cplx> CCx(nx), DDx(nx), ddx(nx);
   std::vector<cplx> CCy(ny), DDy(ny), ddy(ny);
   std::vector<cplx> Uc(N);

   std::vector<double> sigX(nx, 0.0);
   std::vector<double> sigY(ny, 0.0);
   std::vector<cplx>   sX(nx), sY(ny);
   std::vector<cplx>   alpha(nx), alphaP(nx);
   std::vector<cplx>   beta(ny), betaP(ny);

   for(int i=0; i<Lx; ++i) 
      sigX[i]=sig0 * (Lx-i) * (Lx-i) / (1.0*Lx*Lx);
   for(int i=nx-Lx; i<nx; ++i)
      sigX[i]=sig0*(i-nx+Lx)*(i-nx+Lx)/(1.0*Lx*Lx);
   for(int j=0; j<Ly; ++j) 
      sigY[j]=sig0*(Ly-j)*(Ly-j)/(1.0*Ly*Ly);
   for(int j=ny-Ly; j<ny; ++j) 
      sigY[j]=sig0*(j-ny+Ly)*(j-ny+Ly)/(1.0*Ly*Ly);
   for(int i=0; i<nx; ++i) 
      sX[i]=1.0 + I*sigX[i];
   for(int j=0; j<ny; ++j) 
      sY[j]=1.0 + I*sigY[j];
   for(int i=1; i<nx-1; ++i) {
      alpha[i]  = -I*dz/(2.0*ks*dx*dx*sX[i]*(sX[i]+sX[i-1]));
      alphaP[i] = -I*dz/(2.0*ks*dx*dx*sX[i]*(sX[i]+sX[i+1]));
   }
   for(int j=1; j<ny-1; ++j) {
      beta[j]   = -I*dz/(2.0*ks*dy*dy*sY[j]*(sY[j]+sY[j-1]));
      betaP[j]  = -I*dz/(2.0*ks*dy*dy*sY[j]*(sY[j]+sY[j+1]));
   }

   //------------ Field Update ADI -----------------
   for(int h=0; h<numHarmony; ++h)  
   {
      double H = static_cast<double>(D.harmony[h]);
      double invH = 1.0/H;
      
      for(int sliceI=startI; sliceI<endI; ++sliceI) {
         // first step : Y-direction implicit
         // Boundary condition
         for(int j=0; j<ny; ++j) {
            Uc[j*nx + 0]      = 0.0 + I*0.0;
            Uc[j*nx + (nx-1)] = 0.0 + I*0.0;
         }
         for(int i=0; i<nx; ++i) {
            Uc[0*nx + i]        = 0.0 + I*0.0;
            Uc[(ny-1)*nx + i]   = 0.0 + I*0.0;
         }

         for(int j=1; j<ny-1; ++j) {
            // cal. dd
            for(int i=1; i<nx-1; ++i) {
               ddx[i] = - betaP[j] * Un[h][sliceI*N + (j+1)*nx + i]
                        + (H + betaP[j] + beta[j]) * Un[h][sliceI*N + j*nx + i]
                        - beta[j] * Un[h][sliceI*N + (j-1)*nx + i]
                        + Sc[h][sliceI*N + j*nx + i] * currentFlag;
            }
            
            // cal. CC, DD
            CCx[1] = (H-alphaP[1]-alpha[1]) / alphaP[1];
            DDx[1] = ddx[1] / alphaP[1];
            for(int i=2; i<nx-2; ++i) {
               CCx[i] = (H-alphaP[i]-alpha[i]) / alphaP[i]
                       - alpha[i] / (alphaP[i]*CCx[i-1]);
               DDx[i] = ddx[i] / alphaP[i]
                       - alpha[i] * DDx[i-1] / (alphaP[i]*CCx[i-1]);
	         }
            CCx[nx-2] = H-alphaP[nx-2]-alpha[nx-2] - 1.0/CCx[nx-3];
            DDx[nx-2] = ddx[nx-2] - DDx[nx-3]/CCx[nx-3];

            // cal. Uc
	         Uc[j*nx + (nx-2)] = DDx[nx-2] / CCx[nx-2];
	         for(int i=nx-3; i>0; --i) 
	            Uc[j*nx + i] = (DDx[i]-Uc[j*nx + (i+1)]) / CCx[i];
         }

         // second step : X-direction implicit
         for(int i=1; i<nx-1; ++i) {
            // cal. dd
            for(int j=1; j<ny-1; ++j) {
               ddy[j] = - alphaP[i] * Uc[j*nx + (i+1)]
                        + (H + alphaP[i] + alpha[i]) * Uc[j*nx + i]
                        - alpha[i] * Uc[j*nx + (i-1)]
                        + Sc[h][sliceI*N + j*nx + i] * currentFlag;
            }
            // cal. CC, DD
            CCy[1] = (H - betaP[1] - beta[1]) / betaP[1];
            DDy[1] = ddy[1] / betaP[1];
            for(int j=2; j<ny-1; ++j) {
               CCy[j] = (H - betaP[j] - beta[j]) / betaP[j]
                       - beta[j] / (betaP[j]*CCy[j-1]);
               DDy[j] = ddy[j] / betaP[j]
                       - beta[j] * DDy[j-1] / (betaP[j]*CCy[j-1]);
	         }
            CCy[ny-2] = H-betaP[ny-2]-beta[ny-2] - 1.0/CCy[ny-3];
            DDy[ny-2] = ddy[ny-2] - DDy[ny-3]/CCy[ny-3];

            // cal. Un
	         Un[h][sliceI*N + (ny-2)*nx + i] = DDy[ny-2] / CCy[ny-2];
	         for(int j=ny-3; j>0; --j) 
	            Un[h][sliceI*N + j*nx + i] = (DDy[j] - Un[h][sliceI*N + (j+1)*nx + i]) / CCy[j];
         }

      }       //End of for(sliceI)
   }
}

void solve_Sc_3D(Domain &D,int iteration)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int startI = 1; 
   int endI = D.subSliceN+1;
   int numHarmony = D.numHarmony;

   double minX=D.minX;	
   double minY=D.minY;
   int nx=D.nx;   
   int ny=D.ny;
   double dx=D.dx;
   double dy=D.dy;
   double dz=D.dz;
   
   double K0=D.K0;
   double ks=D.ks; 
   double ku=D.ku;
   double dBessel = D.dBessel;
   double ue=D.ue;
   cplx fx = 0.0 +I*0.0;                                                                     cplx fy = 0.0 +I*0.0;
   double w[2]={0.0,0.0}, wx[2]={0.0,0.0}, wy[2]={0.0,0.0};



   std::complex<double> U=(1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);

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

   int N=nx*ny;   
   for(int h=0; h<numHarmony; ++h)
      for(int sliceI=0; sliceI<=endI; ++sliceI)
         for(int j=0; j<N; ++j) {
            D.ScUx[h][sliceI*N + j]=0.0+I*0.0;
            D.ScUy[h][sliceI*N + j]=0.0+I*0.0;
         }

   double coef=dz*eCharge*eCharge*mu0*K0*std::sqrt(1.0+ue*ue)/(4.0*ks*eMass*(D.lambda0*D.numSlice)*dx*dy);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);

   int s=0;
   for(auto& LL : D.loadList) {
      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double macro=p->weight;
            double x=p->x[n]; 
            double y=p->y[n];   
            double theta=p->theta[n]; 
            double gam=p->gamma[n];
            double px=p->px[n];       
            double py=p->py[n];
            double xi2 = std::sqrt(2.0) * K0 / (1.0+K0*K0*0.5*(1.0+ue*ue)) 
                       * std::sqrt((1+K0_alpha)*(px*px+ue*ue*py*py)+(1-K0_alpha)*(px*px*ue*ue+py*py));
            double B=std::atan( 
               2.0*K0_alpha*ue*(px-py) 
               / ((1.0+K0_alpha)*(px+ue*ue*py)+(1-K0_alpha)*(px*ue*ue+py))
            );
            cplx expP=std::exp(I*(B-Phi*0.5));
            cplx expM=std::conj(expP);

            int idxI=(x-minX)/dx;
            int idxJ=(y-minY)/dy;
            wx[1]=(x-minX)/dx-idxI;   wx[0]=1.0-wx[1];
            wy[1]=(y-minY)/dy-idxJ;   wy[0]=1.0-wy[1];	  
            if(idxI>=0 && idxI<nx && idxJ>=0 && idxJ<ny)  {
               for(int h=0; h<numHarmony; ++h)  {
                  int H = D.harmony[h];
                  double dbH = static_cast<double>(H);
                  if(H%2==1)  {  //odd harmony
                     double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
                     int idx=(dbH*xi)/dBessel;
                     w[1]=(dbH*xi/dBessel)-idx; w[0]=1.0-w[1];
                     int order=(H-1)/2;
                     double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=(H+1)/2;
                     double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     fx=sign*(J1-K0_alpha*J2);
                     fy=sign*(J1+K0_alpha*J2);
                  } else {
                     double sign = (H/2 % 2 ==0) ? 1.0 : -1.0;
                     int idx=(dbH*xi)/dBessel;
                     w[1]=(dbH*xi/dBessel)-idx; w[0]=1.0-w[1];
                     int order=(H-2)/2;
                     double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=H/2;
                     double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     order=(H+2)/2;
                     double J3=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                     fx=sign*xi2*dbH*0.5*(expM*(J1-K0_alpha*J2)+expP*(J2-K0_alpha*J3));
                     fy=sign*xi2*dbH*0.5*(expM*(J1+K0_alpha*J2)+expP*(J2+K0_alpha*J3));
                  }
                  cplx macro_expTheta_coef=dbH*macro*coef*std::exp(-I*dbH*(theta+Phi*0.5))/gam;
                  for(int ii=0; ii<2; ++ii)
                     for(int jj=0; jj<2; ++jj) { 
                        D.ScUx[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]
                            -= etaX * wx[ii] * wy[jj] * fx * macro_expTheta_coef;
                        D.ScUy[h][sliceI*N + (idxJ+jj)*nx + (idxI+ii)]
                            -= etaY * wx[ii] * wy[jj] * fy * macro_expTheta_coef;
                     }
               }		//End of harmony
            }	//End of if(idxI,idxJ)
         }	         //End of for(n)
      }		//End of for(sliceI)

      /*
      //Calculate Ez space charge

      for(i=0; i<=endI; i++)
         for(j=0; j<nr; j++) 
            for(l=0; l<L; l++)
               for(f=0; f<F; f++)
                  D->Ez[i][j][l][f]=0.0+I*0.0;

      if(D->SCONOFF == OFF) ;
      else {
         coef=eCharge*velocityC*velocityC*mu0*ku/(1+K0*K0*0.5)/M_PI/dr/dr/(D->lambda0*D->numSlice);

         for(i=startI; i<endI; i++)
         {
            for(j=0; j<nr; j++)
               for(l=0; l<L; l++)
                  for(f=0; f<F; f++)
                     Sc[j][l][f]=0.0+I*0.0;
		
            p=D->particle[i].head[s]->pt;
            while(p) {
               macro=p->weight;
               for(n=0; n<numInBeamlet; n++) {
                  x=p->x[n];         y=p->y[n];   
                  theta=p->theta[n]; gam=p->gamma[n];
                  if(x==0) phi = 0;
                  else     phi = atan2(y,x);

                  K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                  r = sqrt(x*x+y*y);
                  j=(int)(r/dr+0.5);
                  wy[1]=r/dr-(int)(r/dr);  wy[0]=1.0-wy[1];
                  if(j>0 && j<nr) {
                     for(l=0; l<L; l++) 
                        for(f=0; f<F; f++) {
                           coefComp=I*coef*cexp(-I*(l+1)*theta-I*f*phi)*macro*(1+K*K*0.5)*(l+1)/(2.0*j);
                           for(jj=0; jj<2; jj++) Sc[j+jj][l][f]+=coefComp*wy[jj];
                        }
                  }
               }    //End of for(n)            
               p=p->next;
            }       //End of while(p)

            // recalculate 
            for(l=0; l<L; l++)
               for(f=0; f<F; f++)  {
                  j = nr-1;
                  y = j*dr; x=0.0;			      
                  K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                  alpha = 2.0*(l+1)*(l+1)*k0*ku*(1+K*K*0.5)/(1+K0*K0*0.5);
                  Lc = -1.0/(dr*dr*j)*(2*j + f*f*log((j+0.5)/(j-0.5))) - alpha;
                  Lm = 1.0/(dr*dr*j)*(j-0.5);
                  cc[j]=Lm/Lc;
                  dd[j]=Sc[j][l][f];

                  for(j=nr-2; j>=1; j--)  {
                     y = j*dr; x=0.0;			      
                     K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                     alpha = 2.0*(l+1)*(l+1)*k0*ku*(1+K*K*0.5)/(1+K0*K0*0.5);
                     Lp = 1.0/(dr*dr*j)*(j+0.5);
                     Lc = -1.0/(dr*dr*j)*(2*j + f*f*log((j+0.5)/(j-0.5))) - alpha;
                     cc[j]=Lm/(Lc-Lp*cc[j+1]);
                     dd[j]=(Sc[j][l][f]-Lp*dd[j+1])/(Lc-Lp*cc[j+1]);
                  }

                  j=0;
                     y = j*dr; x=0.0;			      
                     K=K0*(1.0+ku*ku*0.5*(x*x+y*y));
                     alpha = 2.0*(l+1)*(l+1)*k0*ku*(1+K*K*0.5)/(1+K0*K0*0.5);
                     Lp = 2.0/(dr*dr);
                     Lc = -2.0/(dr*dr) - alpha;
                     cc[j]=0.0;
                     dd[j]=(Sc[j][l][f]-Lp*dd[j+1])/(Lc-Lp*cc[j+1]);

                  j=0;
                     D->Ez[i][j][l][f] = dd[j];
                  for(j=1; j<nr; j++)
                     D->Ez[i][j][l][f] = dd[j]-cc[j]*D->Ez[i][j-1][l][f];
	       }   //End of for(f)
         }		//End of for(i)

      }   //End of if(SCONOFF==ON)
      */  
      s++;
   }

}





void solve_Field_U_1D(Domain &D,int iteration)
{
   int startI = 1;
   int endI = D.subSliceN+1;

   // field update
   for(int h=0; h<D.numHarmony; ++h)
      for(int sliceI=startI; sliceI<endI; ++sliceI) {
         D.Ux[h][sliceI] += D.ScUx[h][sliceI];
         D.Uy[h][sliceI] += D.ScUy[h][sliceI];
      }
}

void solve_Sc_1D(Domain &D,int iteration)
{
   int myrank, nTasks,rank;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   double ue=D.ue;
   cplx U = (1.0-ue*ue + I*2.0*ue)/(1.0+ue*ue);
   double Phi = std::arg(U);

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

   int startI=1; 
   int endI=D.subSliceN+1;
   double dz=D.dz;
   double ks=D.ks; 
   double ku=D.ku;
   int numHarmony=D.numHarmony;
   double K0=D.K0;
   double dBessel=D.dBessel;
   cplx fx = 0.0 +I*0.0;
   cplx fy = 0.0 +I*0.0;
   double w[2]={0.0,0.0};
   //double xi = ks/ku*K0*K0*(1.0-ue*ue)/(8.0*gam*gam);
   double xi = 0.25*K0*K0*(1.0-ue*ue)/(1.0+(1.0+ue*ue)*K0*K0*0.5);

   int s=0;
   for(auto& LL : D.loadList) {
      double emitX=LL.emitX/D.gamma0;
      double emitY=LL.emitY/D.gamma0;
      double gammaX=(1+LL.alphaX*LL.alphaX)/LL.betaX;
      double gammaY=(1+LL.alphaY*LL.alphaY)/LL.betaY;   
      double sigX = std::sqrt(emitX/gammaX);
      double sigY = std::sqrt(emitY/gammaY);   
      double area = 2.0*M_PI*sigX*sigY;
     
      double coef = dz*eCharge*eCharge*mu0*K0*std::sqrt(1.0+ue*ue)/(2.0*ks*eMass*(D.lambda0*D.numSlice)*area);

      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double th=p->theta[n];      
            double macro=p->weight;
            double gam = p->gamma[n]; 

            for(int h=0; h<numHarmony; h++)  {
               //if (h >= static_cast<int>(D.ScUx.size())) continue;

               int H = D.harmony[h];
               double sign = ((H-1)/2 % 2 ==0) ? 1.0 : -1.0;
               if(H%2==1)  {  //odd harmony
                  int idx=(int)(H*xi/dBessel);
                  w[1]=(H*xi/dBessel)-idx; w[0]=1.0-w[1];
                  int order=(H-1)/2;
                  double J1=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                  order=(H+1)/2;
                  double J2=D.BesselJ[idx][order]*w[0]+D.BesselJ[idx+1][order]*w[1];
                  fx=sign*(J1-K0_alpha*J2);
                  fy=sign*(J1+K0_alpha*J2);
               } else {
                  fx=0.0+I*0.0;
                  fy=0.0+I*0.0;
               }

               cplx macro_expTh_coef=macro*coef*std::exp((-1.0*H)*I*(th+Phi*0.5))/gam;
               D.ScUx[h][sliceI]-=etaX*fx*macro_expTh_coef;
               D.ScUy[h][sliceI]-=etaY*fy*macro_expTh_coef;
            }
         }  //End of for(particle cnt)
         
      }
      
      s++;
   }

   // Calculate Ez space charge
	int L = D.SCLmode;

   for(int l=0; l<L; ++l)
      for(int sliceI=0; sliceI<=endI; ++sliceI)
         D.Ez[l][sliceI]=cplx{0.0,0.0};

   s = 0;
   for (const auto& LL : D.loadList)  
   {
	   int numInBeamlet=LL.numInBeamlet;
      double emitX=LL.emitX/D.gamma0;
      double emitY=LL.emitY/D.gamma0;
      double gammaX=(1+LL.alphaX*LL.alphaX)/LL.betaX;
      double gammaY=(1+LL.alphaY*LL.alphaY)/LL.betaY;   
      double sigX=std::sqrt(emitX/gammaX);
      double sigY=std::sqrt(emitY/gammaY);
   
      double area=2.0*M_PI*sigX*sigY;
      double coef=eCharge*velocityC*velocityC*mu0/ks/area/(D.lambda0*D.numSlice);		 

      for(int sliceI=startI; sliceI<endI; ++sliceI)
      {
         auto& p = D.particle[sliceI].head[s]->pt;
         const size_t cnt=p->x.size();
         for(size_t n=0; n<cnt; ++n) {
            double th=p->theta[n];
            double macro=p->weight;
             
            for(int l=0; l<L; ++l) 
               D.Ez[l][sliceI]+=-1.0*I*coef/(1.0+l)*macro*std::exp(-I*(l+1.0)*th);
         }        
	   }		// End of for(sliceI)
      
      s++;
   }

}

/*    Matrix method but it is very slow.
void solve_Field_U_3D(Domain *D,int iteration)
{
   int h,H,numHarmony,i,j,sliceI,startI,endI,ii;  
   int n,nx,ny;
   double ks,dx,dy,dz,currentFlag;
   double complex alpha,beta,invR,diagB,compVal,*rList,**B,*SList;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   startI=1;   endI=D->subSliceN+1;
   numHarmony=D->numHarmony;

   dx=D->dx;  dy=D->dy;  dz=D->dz;
   ks=D->ks;
   nx=D->nx;  ny=D->ny;
	currentFlag=D->currentFlag;

   // first step
   rList=(double complex *)malloc(nx*sizeof(double complex));
   SList=(double complex *)malloc(nx*sizeof(double complex));
   B=(double complex **)malloc(nx*sizeof(double complex *));
   for(i=0; i<nx; i++)
      B[i]=(double complex *)malloc(nx*sizeof(double complex));

   for(h=0; h<numHarmony; h++)  {
      H = D->harmony[h];
      alpha=-I*dz*0.25/(H*ks)/dx/dx;
  	   beta=-I*dz*0.25/(H*ks)/dy/dy;

      diagB = (1-2*alpha)/alpha;
	   rList[0]=1;
	   rList[1]=-diagB;
      for(i=2; i<nx; i++)
	      rList[i] = -1*(diagB*rList[i-1]+rList[i-2]);
      invR = 1.0/(diagB*rList[nx-1]+rList[nx-2]);

      for(i=0; i<nx; i++)
         for(j=i; j<nx; j++) {
	         compVal = rList[i]*rList[nx-1-j];
            B[i][j] = compVal*invR;
            B[j][i] = compVal*invR;
		   }
	
      for(sliceI=startI; sliceI<endI; sliceI++) {
//       memcpy(&(D->tmpU[0]),&(D->U[h][sliceI][0]),nx*ny*sizeof(double complex ));
         j=0;
            for(i=0; i<nx; i++)
               SList[i]=((1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*D->U[h][sliceI][(j+1)*nx+i]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/alpha;
            for(i=0; i<nx; i++) {
		         compVal=0+I*0;
               for(ii=0; ii<nx; ii++) compVal+=B[i][j]*SList[ii];
               D->Uc[h][sliceI][j*nx+i]=compVal;
            }
		 
         for(j=1; j<ny-1; j++) {
            for(i=0; i<nx; i++)
               SList[i]=((1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*(D->U[h][sliceI][(j-1)*nx+i]+D->U[h][sliceI][(j+1)*nx+i])+D->ScU[h][sliceI][j*nx+i]*currentFlag)/alpha;
            for(i=0; i<nx; i++) {
  	            compVal=0+I*0;
               for(ii=0; ii<nx; ii++) compVal+=B[i][j]*SList[ii];
               D->Uc[h][sliceI][j*nx+i]=compVal;
		      }
         }

         j=ny-1;
            for(i=0; i<nx; i++)
               SList[i]=((1+2*beta)*D->U[h][sliceI][j*nx+i]-beta*D->U[h][sliceI][(j-1)*nx+i]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/alpha;
            for(i=0; i<nx; i++) {
  	            compVal=0+I*0;
               for(ii=0; ii<nx; ii++) compVal+=B[i][j]*SList[ii];
               D->Uc[h][sliceI][j*nx+i]=compVal;
		      }
	   }
   }
   free(rList);
   free(SList);
   for(i=0; i<nx; i++) free(B[i]);
	free(B);

   // second step
   rList=(double complex *)malloc(ny*sizeof(double complex));
   SList=(double complex *)malloc(ny*sizeof(double complex));
   B=(double complex **)malloc(ny*sizeof(double complex *));
   for(i=0; i<ny; i++)
      B[i]=(double complex *)malloc(ny*sizeof(double complex));

   for(h=0; h<numHarmony; h++)  {
	   H = D->harmony[h];
      alpha=-I*dz*0.25/(H*ks)/dx/dx;
  	   beta=-I*dz*0.25/(H*ks)/dy/dy;

      diagB = -1*(1+2*beta)/beta;
	   rList[0]=1;
	   rList[1]=-diagB;
      for(i=2; i<ny; i++)
	      rList[i] = -1*(diagB*rList[i-1]+rList[i-2]);
      invR = 1.0/(diagB*rList[ny-1]+rList[ny-2]);

      for(i=0; i<ny; i++)
         for(j=i; j<ny; j++) {
		      compVal = rList[i]*rList[ny-1-j];
			   B[i][j] = compVal*invR;
			   B[j][i] = compVal*invR;
		   }
	
      for(sliceI=startI; sliceI<endI; sliceI++) {
//       memcpy(&(D->tmpU[0]),&(D->U[h][sliceI][0]),nx*ny*sizeof(double complex ));
         i=0;
            for(j=0; j<ny; j++)
               SList[j]=((1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*D->Uc[h][sliceI][j*nx+(i+1)]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/beta;
            for(j=0; j<ny; j++) {
		         compVal=0+I*0;
               for(ii=0; ii<ny; ii++) compVal+=B[i][j]*SList[ii];
               D->U[h][sliceI][j*nx+i]=compVal;
		      }		 
         for(i=1; i<nx-1; i++) {
            for(j=0; j<ny; j++)
               SList[j]=((1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*(D->Uc[h][sliceI][j*nx+(i-1)]+D->Uc[h][sliceI][j*nx+(i+1)])+D->ScU[h][sliceI][j*nx+i]*currentFlag)/beta;
            for(j=0; j<ny; j++) {
  	            compVal=0+I*0;
               for(ii=0; ii<ny; ii++) compVal+=B[i][j]*SList[ii];
               D->U[h][sliceI][j*nx+i]=compVal;
		      }
         }
         i=nx-1;
            for(j=0; j<ny; j++)
               SList[j]=((1+2*alpha)*D->Uc[h][sliceI][j*nx+i]-alpha*D->Uc[h][sliceI][j*nx+(i-1)]+D->ScU[h][sliceI][j*nx+i]*currentFlag)/beta;
            for(j=0; j<ny; j++) {
  	            compVal=0+I*0;
					for(ii=0; ii<ny; ii++) compVal+=B[i][j]*SList[ii];
               D->U[h][sliceI][j*nx+i]=compVal;
		      }
	   }
   }

   free(rList);
   free(SList);
   for(i=0; i<ny; i++) free(B[i]);
	free(B);	
}
*/
