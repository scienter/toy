#include "particle.h"
#include <complex>

using cplx = std::complex<double>;

enum class OperationMode {
    Unknown       = 0,
    Static        = 1,
    Time_Dependent = 2,
    Twiss         = 3,
    // 필요하면 다른 모드 추가
};



typedef struct _Domain
{
   int dimension;
   double ks,lambda0;
   OperationMode mode;
   int sliceN;   

   int numHarmony, *harmony;
}  Domain; 


typedef struct _Particle 
{
   // Particle List Header
   ptclHead **head;            
}  Particle;


typedef struct _UndulatorList  {
   int undType,numbers,air;
   double *unitStart,*unitEnd;
   double *undStart,*undEnd;
   double *K0,ue;
   int alpha;     //K0_alpha 1:By 0:Bx
   double lambdaU;
   double taper,linTaper,quadTaper;

   struct _UndulatorList *next;
} UndulatorList;
   
typedef struct _QuadList  {
   int numbers;
   double *unitStart,*unitEnd;
   double *qdStart,*qdEnd;
   double *g;	//[T/m]

   struct _QuadList *next;
} QuadList;

typedef struct _PhaseShifter  {
   int num, *step;
   double phase;

   struct _PhaseShifter *next;
} PhaseShifter;


typedef struct _ChiList  {
   int chiON;
	double chiStart,chiEnd,ld,L1,L2,B0,delay,shiftY;

	//self seeding
	int selfSeedON,noiseONOFF,washONOFF,type;
	double d, bragTh, extincL, rangeE,shiftE;
	cplx chi0;
	
	struct _ChiList *next;
} ChiList;


void parameterSetting(Domain *D,const char *input);
int FindParameters (const char *block, int rank, const char *options, const char *input, char *ret);
