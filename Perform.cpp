#ifdef PATHSIM_TESTMODE

#include "Perform.h"

namespace PathSim {

#define K_2PI ( 8.0 * atan(1.0) )		// 2 Pi
#define K_HIST_RES 500	

static double sum = 0.0;
static double rms = 0.0;
static int cnt = 0;
static bool delay = false;
static FILE *stream = NULL;
static double timeinc = 0.0;
static int HistArray[K_HIST_RES];
static double endfreq;
static double freqinc;

static double testfreq = -1.0;

double gDebug1 = 0.0;
double gDebug2 = 0.0;
int iDebug3 = 0;

static LONGLONG StartTime;
static LONGLONG StopTime;
static LONGLONG DeltaTime;
static LONGLONG CountFreq;
static LONGLONG DeltaTimeMax;
static LONGLONG DeltaTimeMin;
static LONGLONG DeltaTimeAve;
static LONGLONG DeltaSamples;

// call to initialize the prformance timer
void InitPerformance()
{
	QueryPerformanceFrequency( (LARGE_INTEGER*)&CountFreq );	//get clock freq
	DeltaTimeMax = 0;
	DeltaTimeAve = 0;
	DeltaSamples = 0;
	DeltaTimeMin = 0x7FFFFFFFFFFFFFFF;
}

// Starts the performance timer
void StartPerformance()
{
	QueryPerformanceCounter( (LARGE_INTEGER*)&StartTime );
}

// Stop performance timer and calculate timing values
void StopPerformance()
{
	QueryPerformanceCounter( (LARGE_INTEGER*)&StopTime );
	DeltaTime = StopTime-StartTime;
	DeltaTimeAve += DeltaTime;
	DeltaSamples++;
	if( DeltaTime>DeltaTimeMax )
		DeltaTimeMax = DeltaTime;
	if( DeltaTime<DeltaTimeMin )
		DeltaTimeMin = DeltaTime;

}

// Call this to measure time between succesive calls to SamplePerformance()
void SamplePerformance()
{
	if(	DeltaSamples == 0 )
	{
		QueryPerformanceCounter( (LARGE_INTEGER*)&StartTime );
	}
	else
	{
		QueryPerformanceCounter( (LARGE_INTEGER*)&StopTime );
		DeltaTime = StopTime-StartTime;
		DeltaTimeAve += DeltaTime;
		if( DeltaTime>DeltaTimeMax )
			DeltaTimeMax = DeltaTime;
		if( DeltaTime<DeltaTimeMin )
			DeltaTimeMin = DeltaTime;
		QueryPerformanceCounter( (LARGE_INTEGER*)&StartTime );
	}
	DeltaSamples++;
}

// output various timing statistics to Message Box
// DeltaTimeMax == maximum time between start()-stop() or sample()-Sample()
// DeltaTimeMin == minimum time between start()-stop() or sample()-Sample()
// DeltaTimeAve == average time between start()-stop() or sample()-Sample()
// DeltaSamples == number of time samples captured
void ReadPerformance()
{
char buf[200];
	if(DeltaSamples != 0 )
	{
		DeltaTime = (DeltaTime*1000000)/CountFreq;
		DeltaTimeMin = (DeltaTimeMin*1000000)/CountFreq;
		DeltaTimeMax = (DeltaTimeMax*1000000)/CountFreq;
		DeltaTimeAve = DeltaTimeAve/DeltaSamples;
		DeltaTimeAve = (DeltaTimeAve*1000000)/CountFreq;
		sprintf( buf, " Max=%I64u uSec Min=%I64u uSec\nAve=%I64u uSec #Samps=%I64u",
					DeltaTimeMax,DeltaTimeMin,DeltaTimeAve,DeltaSamples);
		AfxMessageBox( buf );
	}
}

/////////////////////////////////////////////////////////////////////////
//////  P a t h S i m   T e s t   B e n c h  R o u t i n e s ////////////
/////////////////////////////////////////////////////////////////////////

// Measure RMS power of complex sample over 'samplength' samples
//  from sweep generator.  Save results into a file.
void CalcCpxSweepRMS(cmplx sample, int samplength)
{
	if(testfreq<0.0)	//if first call then ignor this call
		return;
	if( delay )		//allow things to settle before starting sampling
	{
		if(++cnt >= samplength )
		{
			delay = false;
			cnt = 0;
		}
	}
	else
	{
		sum = sum + ( sample.x*sample.x + sample.y*sample.y);
		if(++cnt >= samplength)
		{
			rms = sqrt(sum/(double)samplength);
			sum = 0.0;
			cnt = 0;
			delay = true;		// delay one run to allow things to settle
			if(testfreq <= endfreq)
			{
				fprintf( stream, "%g, %g\n", testfreq, 20*log10(rms + 1e-100) );
				testfreq += freqinc;
				gDebug1 = rms;
				gDebug2 = testfreq;
			}
		}
	}
}
// Measure RMS power of real sample over 'samplength' samples
//  from sweep generator.  Save results into a file.
void CalcSweepRMS(double sample, int samplength)
{
	if(testfreq<0.0)	//if first call then ignor this call
		return;
	if( delay )		//allow things to settle before starting sampling
	{
		if(++cnt >= samplength )
		{
			delay = false;
			cnt = 0;
		}
	}
	else
	{
		sum = sum + (sample*sample);
		if(++cnt >= samplength)
		{
			rms = sqrt(sum/(double)samplength);
			sum = 0.0;
			cnt = 0;
			delay = true;		// delay one run to allow things to settle
			if(testfreq <= endfreq)
			{
				fprintf( stream, "%g, %g\n", testfreq, 20*log10(rms) );
				testfreq += freqinc;
				gDebug1 = rms;
				gDebug2 = testfreq;
			}
		}
	}
}

// Measure phase angle of complex sinwave sample over 'samplength' samples
//  from sweep generator.  Save results into a file.
// Uses trig identity and assumes sinwave inputs of same frequency:
//  sin(a)sin(b) = 1/2 cos(a-b) - 1/2 cos(a+b)
// LP filter(average) the cos(a+b) term out and it leaves 1/2cos(phzedif)
//    then take inv cos to get angle.
void CalcCpxSweepPhz(cmplx sample, int samplength)
{
double mag;
	if(testfreq<0.0)	//if first call then ignor this call
		return;
	if( delay )		//allow things to settle before starting sampling
	{
		if(++cnt >= samplength )
		{
			delay = false;
			cnt = 0;
		}
	}
	else
	{	//normalize amplitude.
		mag = (sample.x*sample.x+sample.y*sample.y);
		sum = sum + sample.x*sample.y/mag;	//add sin(a)*sin(b) term
		if(++cnt >= samplength)
		{
			rms = 360.0*acos(2*sum/(double)samplength)/K_2PI; //convert to angle
			sum = 0.0;
			cnt = 0;
			delay = true;		// delay one run to allow things to settle
			if(testfreq <= endfreq)
			{
				fprintf( stream, "%g, %g\n", testfreq, rms);
				testfreq += freqinc;
				gDebug1 = rms;
				gDebug2 = testfreq;
			}
		}
	}
}

// Generate a RMS=1 real sinwave sweep frequency from start to stop
//   with given sample rate and step size.
void SweepGen(  double* output, double samprate, 
							double start, double stop, double step )
{
	if(testfreq<0.0)	//if first call then initialize everything
	{
		timeinc = 0.0;
		testfreq = start;
		endfreq = stop;
		freqinc = step;
		sum = 0.0;
		rms = 0.0;
		cnt = 0;
		delay = true;		// delay one run to allow things to settle
		stream = fopen( "c:\\Data.prn", "wt" );	//Open file for write
	}
	*output = 1.414213562*sin(timeinc);
	timeinc += ( (K_2PI/samprate)*testfreq);
	timeinc = fmod(timeinc,K_2PI);	//keep radian counter bounded
}
// Generate a RMS=1 complex sinwave sweep frequency from start to stop
//   with given sample rate and step size.
void SweepGenCpx(  cmplx* output, double samprate, 
							double start, double stop, double step )
{
	if(testfreq<0.0)	//if first call then initialize everything
	{
		timeinc = 0.0;
		testfreq = start;
		endfreq = stop;
		freqinc = step;
		sum = 0.0;
		rms = 0.0;
		delay = true;		// delay one run to allow things to settle
		cnt = 0;
		stream = fopen( "c:\\Data.prn", "wt" );	//Open file for write
	}
	output->x = cos(timeinc);
	output->y = sin(timeinc);
	timeinc += ( (K_2PI/samprate)*testfreq);
	timeinc = fmod(timeinc,K_2PI);	//keep radian counter bounded
}


void HistogramSamp(  double sample, double min, double max, int numsamps )
{
double K;
int i;
int hmax;
	if(!delay)	//if first call then initialize everything
	{
		cnt = 0;
		delay = true;		// delay one run to allow things to settle
		stream = fopen( "c:\\Data.prn", "wt" );	//Open file for write
		for(i=0;i<K_HIST_RES; i++)
			HistArray[i] = 0;
	}
	if(cnt < numsamps)
	{
		K = (K_HIST_RES-1)/(max-min);
		i = (int)( K*(sample-min) );
		if( i>0 && i <K_HIST_RES)
			HistArray[i]++;
		cnt++;
		gDebug1 = (double)cnt;
	}
	if(cnt == numsamps)
	{
		hmax = 0;
		for(i=0;i<K_HIST_RES; i++)
		{
			if( HistArray[i] > hmax)
				hmax = HistArray[i];
		}
		for(i=0;i<K_HIST_RES; i++)
		{
			K = (double)i*(max-min)/(double)K_HIST_RES + min;
			fprintf( stream, "%g, %g\n", K,
							(double)HistArray[i]/(double)hmax);
		}
		cnt++; 
	}
	
}

//call this to end the test and close the data file if open;
void EndTest(void)
{
	delay = false;
	timeinc = 0.0;
	testfreq = -1.0;
	if (stream)
		fclose( stream );
}

///////////////////////////////////////////////////////////////////
//calculates RMS value of complex sample over numsamps
//  called once for each sample being measured
///////////////////////////////////////////////////////////////////
double CalcCpxRMS(cmplx sample, int numsamps)
{
static double sum = 0.0;
static double rms = 0.0;
	sum = sum + ( sample.x*sample.x + sample.y*sample.y);
	if(++cnt >= numsamps)
	{
		rms = sqrt(sum/(double)numsamps);
		sum = 0.0;
		cnt = 0;
	}
	return rms;
}

///////////////////////////////////////////////////////////////////
//calculates running average RMS value of buf[] samples
//  called once for each buffer to be measured
///////////////////////////////////////////////////////////////////
double CalcRunAveRMS(double* buf, int bufsize)
{
#define RMSBAVE 1000
static double sum = 0.0;
static double rms = 0.0;
	sum = 0.0;
	for( int i = 0; i<bufsize; i++)
		sum = sum + ( buf[i]*buf[i] );
	rms = (1.0/RMSBAVE)*sqrt(sum/bufsize) + (1.0-1.0/RMSBAVE)*rms;
	return rms;
}

///////////////////////////////////////////////////////////////////
//calculates running average RMS value of buf[] complex samples
//  called once for each buffer to be measured
///////////////////////////////////////////////////////////////////
double CalcCpxRunAveRMS(cmplx* buf, int bufsize)
{
#define RMSCPXAVE 300
static double sum = 0.0;
static double rms = 0.0;
	for( int i = 0; i<bufsize; i++)
		sum = sum + ( buf[i].x*buf[i].x +  buf[i].y*buf[i].y);
	rms = (1.0/RMSCPXAVE)*sqrt(sum/bufsize) + (1.0-1.0/RMSCPXAVE)*rms;
	sum = 0.0;
	return rms;
}

///////////////////////////////////////////////////////////////////
//calculates RMS value of buf[] samples over RMS_SIZE*numsamps
//  called once for each buffer to be measured
///////////////////////////////////////////////////////////////////
double CalcTotalRMS(double* buf, int numsamps)
{
#define RMS_SIZE 7000	//abt 30 minutes at 8000 KSP
static double sum = 0.0;
static int cnt = 0;
static double rms = 0.0;
	for( int i = 0; i<numsamps; i++)
		sum = sum + ( buf[i]*buf[i] );
	if( ++cnt >= RMS_SIZE)
	{
		rms = sqrt(sum/(cnt*numsamps));
		sum = 0.0;
		cnt = 0;
	}
	return rms;
}

} // namespace PathSim

#endif // PATHSIM_TESTMODE
