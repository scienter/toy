#ifndef MESH_H
#define MESH_H

#include <complex>
#include <vector>
#include "particle.h"

using cplx = std::complex<double>;
const std::complex<double> I(0.0,1.0);

enum class ONOFFMode {
    OFF    = 0,
    ON     = 1
};

enum class OperationMode {
    Unknown    = 0,
    Static     = 1,
    Time_Dependent = 2,
    Twiss      = 3
};

enum class UndMode {
    Unknown    = 0,
    BiPolar    = 1,
    QuadPolar   = 2
    // 필요하면 다른 모드 추가
};

struct UndList  {
   UndMode type = UndMode::Unknown;

   double lambdaU=0.0;
   double linTaper,quadTaper,slopeK,ue;
   int K0_alpha=1;    // 1:By 0:Bx
   int numbers=0;
   bool air;

   std::vector<double> unitStart, unitEnd;
   std::vector<double> undStart, undEnd;
   std::vector<double> K0;
};

struct QuadList  {
   int numbers;
   std::vector<double> unitStart, unitEnd;
   std::vector<double> qdStart, qdEnd;
   std::vector<double> g;  //[T/m]
};


typedef struct _Domain
{
   int dimension;
   double ks,lambda0;
   OperationMode mode;

   int numHarmony, *harmony;
   std::vector<std::vector<cplx>> Ux,Uy,ScUx,ScUy;
   std::vector<std::vector<cplx>> Ez;
   std::vector<std::vector<double>> totalEnergyX,totalEnergyY;


   //Save option
   int maxStep;

   //Domain box
   int sliceN,subSliceN,nx,ny,minI,maxI;
   double minX,maxX,minY,maxY,minZ,maxZ;
   double Lz,dz;   
   double dx,dy;

   // ABC condition
   int abcN;
   double abcSig;

   //Electron beam
   int nSpecies=0;
   int numSlice;
   double gamma0=1.0;
   std::vector<LoadList> loadList;
   std::vector<Particle> particle;

   //Undulator
   UndMode undType;
   bool currentFlag,driftFlag;
   int numLambdaU,nUnd;
   double lambdaU,ku,K0=0.0,K0_alpha,ue;
   std::vector<UndList> undList;

   //Quad
   int nQuad;
   double g;
   std::vector<QuadList> quadList;
   
   //Space charge
   int SCLmode,SCFmode;
   double dr;

   //Bessel Table
   double BesselMax,dBessel;
   int BesselMaxOrder,bn;
   std::vector<std::vector<double>> BesselJ;

   //Seed
   double P0,duration,spotSigR,a0,zR,focus;

}  Domain; 




void parameterSetting(Domain *D,const char *input);
int FindParameters (const char *block, int rank, const char *options, const char *input, char *ret);
void boundary(Domain *D);
void loadBeam(Domain *D,LoadList &LL,int s,int iteration);
void updateK_quadG(Domain *D,int iteration,double half);
void push_theta_gamma(Domain *D,int iteration);
void solveField(Domain *D,int iteration);
void updateTotalEnergy(Domain *D,int iteration);
void transversePush(Domain *D,int iteration);
void calculate_twiss(Domain &D,int iteration);

void saveFieldsToTxt(const Domain &D, const std::string& fileName);
void saveParticlesToTxt(const Domain &D, int species, const std::string& fileName);


#endif   //MESH_H
