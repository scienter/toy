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

typedef struct _Domain
{
   int dimension;
   double ks,lambda0;
   OperationMode mode;

   int numHarmony, *harmony;
   std::vector<cplx> Ux,Uy,ScUx,ScUy;

   //Domain box
   int sliceN,subSliceN,nx,ny,minI,maxI;
   double minZ, maxZ;   

   //Electron beam
   int nSpecies=0;
   int numSlice;
   std::vector<LoadList> loadList;
   std::vector<Particle> particle;

   //Undulator
   int numLambdaU;

}  Domain; 

void parameterSetting(Domain *D,const char *input);
int FindParameters (const char *block, int rank, const char *options, const char *input, char *ret);
void boundary(Domain *D);
void loadBeam(Domain *D,LoadList &LL,int s,int iteration);


void saveParticlesToTxt(const Domain &D, int species, const std::string& fileName);

#endif   //MESH_H
