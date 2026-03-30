#ifndef MESH_H
#define MESH_H

#include <complex>
#include <vector>
#include "particle.h"

using cplx = std::complex<double>;

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
    Normal    = 1,
    AppleX   = 2
    // 필요하면 다른 모드 추가
};

struct UndList  {
   UndMode type = UndMode::Unknown;

   double lambdaU=0.0;
   double linTaper,quadTaper,slopeK,ue=0.0;
   int K0_alpha=1;    // 1:By 0:Bx
   int numbers=0;
   bool air;

   std::vector<double> unitStart, unitEnd;
   std::vector<double> undStart, undEnd;
   std::vector<double> K0;


};


typedef struct _Domain
{
   int dimension;
   double ks,lambda0;
   OperationMode mode;

   int numHarmony, *harmony;
   std::vector<cplx> Ux,Uy,ScUx,ScUy;


   //Save option
   int maxStep;

   //Domain box
   int sliceN,subSliceN,nx,ny,minI,maxI;
   double minZ,maxZ;
   double Lz,dz;   

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

}  Domain; 




void parameterSetting(Domain *D,const char *input);
int FindParameters (const char *block, int rank, const char *options, const char *input, char *ret);
void boundary(Domain *D);
void loadBeam(Domain *D,LoadList &LL,int s,int iteration);
void updateK_quadG(Domain *D,int iteration,double half);
void push_theta_gamma(Domain *D,int iteration);
void saveParticlesToTxt(const Domain &D, int species, const std::string& fileName);

#endif   //MESH_H
