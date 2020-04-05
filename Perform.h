#ifndef PATHSIM_PERFORM_HPP
#define PATHSIM_PERFORM_HPP

#include "cmplx.h"

namespace PathSim {

extern void InitPerformance();
extern void StartPerformance();
extern void StopPerformance();
extern void ReadPerformance();
extern void SamplePerformance();

extern double testfreq;
extern double gDebug1;
extern double gDebug2;
extern int iDebug3;

extern void CalcCpxSweepRMS(cmplx sample, int samplength);
extern void CalcSweepRMS(double sample, int samplength);
extern void CalcCpxSweepPhz(cmplx sample, int samplength);

extern void SweepGen( double* output, double samprate, 
							double start, double stop, double step);
extern void SweepGenCpx( cmplx* output, double samprate, 
							double start, double stop, double step);
extern void HistogramSamp(  double sample, double min, double max, int numsamps );

extern void EndTest(void);
extern double CalcCpxRMS(cmplx sample, int numsamps);
extern double CalcRunAveRMS(double* buf, int bufsize);
extern double CalcCpxRunAveRMS(cmplx* buf, int bufsize);
extern double CalcTotalRMS(double* buf, int numsamps);

} // namespace PathSim

#endif // PATHSIM_PERFORM_HPP
