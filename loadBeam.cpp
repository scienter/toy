#include "mesh.h"
#include "constants.h"
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <mpi.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>

void loadBeam1D(Domain &D,LoadList &LL,int s,int iteration);
void loadBeam3D(Domain &D,LoadList &LL,int s,int iteration);

void loadBeam(Domain *D,LoadList &LL,int s,int iteration)
{
   switch(D->dimension)  {
   case 1:
      loadBeam1D(*D,LL,s,iteration);
      break;
   case 3:
      loadBeam3D(*D,LL,s,iteration);
      break;
   default:
      break;
   }
}

void loadBeam3D(Domain &D,LoadList &LL,int s,int iteration)
{
   int myrank,nTasks,rank;
	MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   int startI=1;       
   int endI=1+D.subSliceN;
   
   int numInBeamlet=LL.numInBeamlet;
   double gamma0=LL.energy/mc2+1;
   double dGam=LL.spread*gamma0;
   double current=LL.peakCurrent;		// peak current in a cell
   double bucketZ=D.lambda0*D.numSlice;	// size of a big slice
   double ptclCnt=numInBeamlet*LL.numBeamlet;	

   // Calculation recommanding quad g*l, beta_min, beta_max
   for(auto& QD : D.quadList) {
      double beta0=0.5*(LL.betaX+LL.betaY);
      double L = 0.5*(QD.unitEnd[0]-QD.unitStart[0]);
      double lquad = QD.qdEnd[0]-QD.qdStart[0];
      double gl=2.0*gamma0*eMass*velocityC/eCharge/beta0*sqrt(2.0/(1.0+std::sqrt(1.0+4.0*L*L/beta0/beta0)));
      double g = gl/lquad;
      double tmp = 2.0*gamma0*eMass*velocityC/(eCharge*gl);
      double min_beta = tmp*(tmp/L-1.0)/sqrt(tmp*tmp/L/L-1.0);
      double max_beta = tmp*(tmp/L+1.0)/sqrt(tmp*tmp/L/L-1.0);
      if (myrank==0) 
         printf("Recommandations : quad g=%g, quad K=%g, cen_beta=%g, min_beta=%g, max_beta=%g\n",g,eCharge/(gamma0*eMass*velocityC)*g,beta0,min_beta,max_beta);
   }

   double dPhi=2.0*M_PI*D.numSlice;
   double div=2.0*M_PI/(1.0*numInBeamlet);
   double macro=current/eCharge/velocityC*bucketZ/ptclCnt;

   size_t totalCnt=LL.numBeamlet*LL.numInBeamlet;
   
   // gsl random generator
   gsl_rng_env_setup();

   const gsl_rng_type * T = gsl_rng_default;
   gsl_rng *ran = gsl_rng_alloc(T);
   
   gsl_qrng *q1 = nullptr;

   //q1 = gsl_qrng_alloc(gsl_qrng_niederreiter_2,3);
   q1=gsl_qrng_alloc(gsl_qrng_sobol,7);	

   unsigned long randskip = 0;   
   if (LL.randONOFF==false) {
      srand(static_cast<unsigned int>(myrank));
      randskip = static_cast<unsigned long>(myrank);
   }
   else {
      srand(static_cast<unsigned int>(time(nullptr)));
      randskip = static_cast<unsigned long>(rand()) % (nTasks * 2UL);
   }
   double v1[7] = {0.0};
   for (unsigned long ii = 0; ii < randskip; ++ii) {
      gsl_qrng_get(q1, v1);
   }
   
   for(int sliceI=startI; sliceI<endI; ++sliceI) {
      //position define     
      double posZ=(sliceI-startI+D.minI)*bucketZ+D.minZ;
      double n0=0.0;
      if(LL.type==BeamMode::Polygon) {
         for(int l=0; l<LL.znodes-1; ++l) {
            if(posZ>=LL.zpoint[l] && posZ<LL.zpoint[l+1])
               n0=(LL.zn[l+1]-LL.zn[l])/(LL.zpoint[l+1]-LL.zpoint[l])*(posZ-LL.zpoint[l])+LL.zn[l];
         }
      }
      
      double emitX=LL.emitX/gamma0;
      double emitY=LL.emitY/gamma0;
      double gammaX=(1+LL.alphaX*LL.alphaX)/LL.betaX;
      double gammaY=(1+LL.alphaY*LL.alphaY)/LL.betaY;   
      double sigX=sqrt(emitX/gammaX);
      double sigY=sqrt(emitY/gammaY);
      double sigXPrime=sqrt(emitX*gammaX);
      double sigYPrime=sqrt(emitY*gammaY);

      double distanceX=std::sqrt(std::fabs((LL.betaX-1.0/gammaX)/gammaX));
      double distanceY=std::sqrt(std::fabs((LL.betaY-1.0/gammaY)/gammaY));
      double vz=std::sqrt(gamma0*gamma0-1.0)/gamma0;	//normalized
      double delTX=distanceX/vz;	//normalized
      double delTY=distanceY/vz;	//normalized
      if(vz==0.0) { delTX=delTY=0.0; }

      unsigned int beamlets=LL.numBeamlet*n0;
      double remacro = (beamlets >0)
                     ? macro*static_cast<double>(LL.numBeamlet)/beamlets*n0
                     : 0.0;
      double eNumbers=remacro*numInBeamlet;
      if(eNumbers<10) eNumbers=10;  

   }
/*

   for(i=startI; i<endI; i++) {
	 }

      index=0;
      for(b=0; b<beamlets; b++)  {
         if(index>=D->numSlice) index=0; else ;
         index=0;
         gsl_qrng_get(q1,v1);

         r1=v1[0];  if(r1==0.0)  r1=1e-9;        r2=v1[1];
         pr1=v1[2]; if(pr1==0.0) pr1=1e-9;       pr2=v1[3];
         th=v1[4];	 gam=v1[5];   if(gam==0.0) gam=1e-9;
	     //th+=randTh;
		  //tmpInt=(int)th;
		  //th-=tmpInt;

         if(LL->transFlat==OFF) {        //Transverse Flat
            coef=sqrt(-2.0*log(r1));
            x=coef*cos(2*M_PI*r2);
            x*=sigX;
            y=coef*sin(2*M_PI*r2);
            y*=sigY;

            coef=sqrt(-2.0*log(pr1));
            xPrime=coef*cos(2*M_PI*pr2);
            xPrime*=sigXPrime;
            yPrime=coef*sin(2*M_PI*pr2);
            yPrime*=sigYPrime;
         }  else  {                       //Transverse Gaussian 
            coef=sqrt(r1);
            x=coef*cos(2*M_PI*r2);
	    x*=sigX;
            y=coef*sin(2*M_PI*r2);
	    y*=sigY;

            //coef=sqrt(pr1);
            //xPrime=coef*cos(2*M_PI*pr2);
            //xPrime*=sigXPrime;
            //yPrime=coef*sin(2*M_PI*pr2);
	    //yPrime*=sigYPrime;
            coef=sqrt(-2.0*log(pr1));
            xPrime=coef*cos(2*M_PI*pr2);
            xPrime*=sigXPrime;
            yPrime=coef*sin(2*M_PI*pr2);
            yPrime*=sigYPrime;
         }

         tmp=sqrt(-2.0*log(gam))*cos(v1[6]*2*M_PI);
	 gam=gamma0+sigGam*tmp*ESn0;
	 //coef=sqrt(-2.0*log(v1[6]));
	 //dGam=coef*(2*gam-1)*sigGam;
	 //gam=gamma0+dGam;
	 theta0=th*dPhi;

         pz=sqrt((gam*gam-1.0)/(1.0+xPrime*xPrime+yPrime*yPrime));
         px=xPrime*pz;
         py=yPrime*pz;

         x-=delTX*px/gam;
         y-=delTY*py/gam;

         New = (ptclList *)malloc(sizeof(ptclList));
         New->next = D->particle[i].head[s]->pt;
         D->particle[i].head[s]->pt = New;

         New->weight=remacro;
	 cnt+=remacro;
         New->index=LL->index;  	//index
         New->core=myrank;  	

         New->x=(double *)malloc(numInBeamlet*sizeof(double ));
         New->y=(double *)malloc(numInBeamlet*sizeof(double ));
         New->px=(double *)malloc(numInBeamlet*sizeof(double ));
         New->py=(double *)malloc(numInBeamlet*sizeof(double ));
         New->theta=(double *)malloc(numInBeamlet*sizeof(double ));
         New->gamma=(double *)malloc(numInBeamlet*sizeof(double ));

         for(n=0; n<numInBeamlet; n++)  {		     
            New->x[n] = x;  New->y[n] = y+YOff+shiftY;
	    aveY+=New->y[n]*remacro;
            //New->px[n] = px;  New->py[n] = py+PyOff*gamma0;
            New->px[n] = px;  New->py[n] = py;
            New->gamma[n]=gam;

            theta=theta0+n*div;
            noise=0.0;
            for(m=1; m<=numInBeamlet/2; m++) {
               sigma=sqrt(2.0/eNumbers/(m*m*1.0));	//harmony is 1.
  	       an=gsl_ran_gaussian(ran,sigma);
               bn=gsl_ran_gaussian(ran,sigma);
               noise += an*cos(m*theta)+bn*sin(m*theta);
            }
            tmp=theta + noise*noiseONOFF;
            New->theta[n]=tmp;	     
         }		// End for (n)
         LL->index+=1;
      }	//End for (b)
   }     //End for (i)
   gsl_qrng_free(q1);
   gsl_rng_free(ran);

printf("myrank=%d, cnt=%d\n",myrank,cnt);

   for(rank=1; rank<nTasks; rank++) 
      if(myrank==rank) MPI_Send(&cnt,1,MPI_INT,0,myrank,MPI_COMM_WORLD); else ;
   
   if(myrank==0) {
      for(rank=1; rank<nTasks; rank++) {
         MPI_Recv(&recvInt,1,MPI_INT,rank,rank,MPI_COMM_WORLD,&status);
	 cnt+=recvInt;
      }
   } else ;
   for(rank=1; rank<nTasks; rank++) 
      if(myrank==rank) MPI_Send(&aveY,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); else ;
   
   if(myrank==0) {
      for(rank=1; rank<nTasks; rank++) {
         MPI_Recv(&recvDb,1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
	 aveY+=recvDb;
      }
   } else ;
   if(myrank==0) 
      printf("beam index = %d, beam charge = %g [pC], aveY=%g [m],cnt=%d, numInBeamlet=%d\n",s,cnt*1.602e-7*numInBeamlet,aveY/(cnt*1.0*numInBeamlet),cnt,numInBeamlet);
   else ;
*/
}

void loadBeam1D(Domain &D,LoadList &LL,int s,int iteration)
{

   int myrank,nTasks;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   int numInBeamlet=LL.numInBeamlet;
   double gamma0=LL.energy/mc2+1;
   double dGam=LL.spread*gamma0;

   double current=LL.peakCurrent;		// peak current in a cell
   double bucketZ=D.lambda0*D.numSlice;		     	// size of a big slice
   int ptclCnt=numInBeamlet*LL.numBeamlet;	

   double dPhi=2.0*M_PI*D.numSlice;
   double div=2.0*M_PI/(1.0*numInBeamlet);
   bool noiseONOFF=LL.noiseONOFF;

   double macro=current/eCharge/velocityC*bucketZ/ptclCnt;

   // gsl random generator
   gsl_rng_env_setup();

   const gsl_rng_type * T = gsl_rng_default;
   gsl_rng *ran = gsl_rng_alloc(T);
   
   gsl_qrng *q1 = nullptr;

   //q1 = gsl_qrng_alloc(gsl_qrng_niederreiter_2,3);
   q1=gsl_qrng_alloc(gsl_qrng_sobol,3);	

   unsigned long randskip = 0;   
   if (LL.randONOFF==false) {
      srand(static_cast<unsigned int>(myrank));
      randskip = static_cast<unsigned long>(myrank);
   }
   else {
      srand(static_cast<unsigned int>(time(nullptr)));
      randskip = static_cast<unsigned long>(rand()) % (nTasks * 2UL);
   }
   double v1[3] = {0.0};
   for (unsigned long ii = 0; ii < randskip; ++ii) {
      gsl_qrng_get(q1, v1);
   }

   int minI=D.minI;
   int maxI=D.maxI;
   int startI=1;	   
   int endI=1+D.subSliceN;

   for(int i=startI; i<endI; ++i) {
      //position define     
      double posZ=(i-startI+minI)*bucketZ+D.minZ;
      double n0=0.0;
      if(LL.type==BeamMode::Polygon) {
         for(int l=0; l<LL.znodes-1; ++l) {
            if(posZ>=LL.zpoint[l] && posZ<LL.zpoint[l+1])
               n0=(LL.zn[l+1]-LL.zn[l])/(LL.zpoint[l+1]-LL.zpoint[l])*(posZ-LL.zpoint[l])+LL.zn[l];
         }
      }

      unsigned int beamlets=LL.numBeamlet*n0;
      double remacro = (beamlets >0)
                     ? macro*static_cast<double>(LL.numBeamlet)/beamlets*n0
                     : 0.0;
      double eNumbers=remacro*numInBeamlet;
      if(eNumbers<10) eNumbers=10;  

      size_t totalParticles = beamlets * numInBeamlet;
      auto New = std::make_unique<ptclList>();

      // head[s] is nullptr, the generate new.
      if (D.particle[i].head[s] == nullptr) {
         D.particle[i].head[s] = new ptclHead{};
         D.particle[i].head[s]->pt = nullptr;
      }

      New->next = D.particle[i].head[s]->pt;
      D.particle[i].head[s]->pt = New.get();

      New->weight = remacro;
      New->x.resize(totalParticles);
      New->y.resize(totalParticles);
      New->px.resize(totalParticles);
      New->py.resize(totalParticles);
      New->theta.resize(totalParticles);
      New->gamma.resize(totalParticles);
      New->index.resize(totalParticles);
      New->core.resize(totalParticles);

      unsigned long ptclIdx=0;
      for(unsigned int b=0; b<beamlets; ++b)  {
         gsl_qrng_get(q1,v1);
         double th=v1[0];           
         double gam= (v1[1]==0.0) ? 1e-4 : v1[1];
         double tmp=std::sqrt(-2.0*std::log(gam))*std::cos(v1[2]*2.0*M_PI);
         gam=gamma0+dGam*tmp;
         double theta0=th*(dPhi-(numInBeamlet-1.0)/numInBeamlet*1.0 * 2*M_PI);
         //theta0=(th)*(2*M_PI-(numInBeamlet-1.0)/(numInBeamlet*1.0)*2*M_PI);
         //theta0=(th)*(dPhi);  
         for(int n=0; n<numInBeamlet; ++n)  { 
            unsigned long idx = ptclIdx++;
      
            New->x[idx]=0.0;
            New->y[idx]=0.0;
            New->px[idx]=0.0;        
            New->py[idx]=0.0;
            New->gamma[idx]=gam;		//gamma
         
            double theta=theta0+n*div;
            double noise=0.0;
            for(int m=1; m<=numInBeamlet/2; ++m) {
               double sigma=std::sqrt(2.0/eNumbers/(1.0*m*m)); //Fawley PRSTAB V5 070701 (2002)
               //an=gaussianDist_1D(sigma);
               //bn=gaussianDist_1D(sigma);
  	            double an=gsl_ran_gaussian(ran,sigma);
               double bn=gsl_ran_gaussian(ran,sigma);
               noise += an*std::cos(m*theta) + bn*std::sin(m*theta);
            }				 
            double final_theta = theta + noise * (noiseONOFF ? 1.0 : 0.0);
            //if(tmp>=dPhi) tmp-=dPhi; 
            //else if(tmp<0) tmp+=dPhi; 
            //else ;
            New->theta[idx]=final_theta;       
         }  	// End for(n)
      }      // End for(b) 
      New.release();
   }			//End of for(i)

   gsl_qrng_free(q1);
   gsl_rng_free(ran);
}

/*
void random_2D(double *x,double *y,gsl_qrng *q1)
{
   double v[2];

   gsl_qrng_get(q1,v);
   *x=v[0];
   *y=v[1];
}

double gaussianDist_1D(double sigma)
{
   double r,prob,v,z,random;
   int intRand,randRange=1e4;

   r=1.0;
   prob=0.0;
//   gsl_qrng *q=gsl_qrng_alloc(gsl_qrng_niederreiter_2,1);
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
//      gsl_qrng_get(q,&v);
      z = 4.0*(random-0.5);	//up to 3 sigma
      prob=exp(-z*z);
   }
//   gsl_qrng_free(q); 
  
   return z*sigma;
}

double randomValue(double beta)
{
   double r;
   int intRand, randRange=1000, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((double)intRand)/randRange+(1.0-beta);

   return r;
}
*/
