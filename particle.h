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
   double weight;
   int index,core;
   struct ptclList *next = nullptr;
};

struct LoadList  {
   BeamMode type = BeamMode::Unknown;
	int index=0;


   int numInBeamlet;    // For particles in one slice not a big slice.
   int numBeamlet;   // Number of beamlets in one big slice.


   // For Polygon mode
   int znodes = 0;
   std::vector<double> zpoint;
   std::vector<double> zn;

};
   

#endif
