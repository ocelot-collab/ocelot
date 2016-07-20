#ifndef UNDULATORH
#define UNDULATORH
//extern "C" __declspec(dllexport) int trajectory(double *aMagField,
extern "C" int trajectory(double *aMagField,
                                                const int colsMF,
                                                const int lenMF,
                                                const double *misalign,
                                                const int bRough,
                                                const double *aInitCond,
                                                const double *undul_param,
                                                const int NstepMotion,
                                                const int Nsuperperiod,
                                                double *aMotion);
#endif
