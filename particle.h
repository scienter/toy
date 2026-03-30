#ifndef PARTICLE_H
#define PARTICLE_H


#include <vector>

enum class BeamMode {
    Unknown    = 0,
    Polygon    = 1,
    Gaussian   = 2
    // 필요하면 다른 모드 추가
};

struct ptclHead  {
    struct ptclList *pt = nullptr;
};

struct Particle 
{
   std::vector<ptclHead*> head; // species 개수만큼 ptclHead* 를 가짐
};


struct ptclList  {
   std::vector<double> x;
   std::vector<double> y;
   std::vector<double> px;
   std::vector<double> py;
   std::vector<double> theta;
   std::vector<double> gamma;
   std::vector<int> index;
   std::vector<int> core;
   double weight;
   struct ptclList *next = nullptr;
};

struct LoadList  {
   BeamMode type = BeamMode::Unknown;

   int numInBeamlet;    // For particles in one slice not a big slice.
   int numBeamlet;   // Number of beamlets in one big slice.
   bool noiseONOFF, randONOFF;
   
   double energy,spread,gamma0;
   double peakCurrent;

   // For Polygon mode
   int znodes = 0;
   std::vector<double> zpoint;
   std::vector<double> zn;

};
   

#endif
