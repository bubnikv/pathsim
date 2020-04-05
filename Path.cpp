#include "Path.h"

#define K_2PI ( 8.0 * atan(1.0) )			// 2 Pi
#define OFFSET_FREQ_CONST (K_2PI/8000.0)	//2Pi/8000
#define KGNB 0.62665707	//equivalent Noise BW of Gaussian shaped filter

#define RATE_12_8 0		//Used for 0.1 > Spread >= 0.4
#define RATE_64 1		//Used for 0.4 > Spread >= 2.0
#define RATE_320 2		//Used for 2.0 > Spread >= 10.0

namespace PathSim {

Path::Path()
{
    m_Indx = 0;
    m_NoiseSampRate = RATE_320;
}

void Path::InitPath( double Spread, double Offset, int blocksize, int numpaths, bool active)
{
    m_BlockSize = blocksize;
    m_Offset = Offset;
    m_Spread = Spread;
    m_PathActive = active;
    m_FirState0 = INTP_QUE_SIZE-1;
    m_FirState1 = INTP_QUE_SIZE-1;
    m_FirState2 = INTP_QUE_SIZE-1;
    m_FirState3 = INTP_QUE_SIZE-1;
    m_Indx = 0;
    m_inc = 0;
    m_Timeinc = 0.0;
    if( (m_Spread > 2.0) && (m_Spread <= 30.0) )
    {
        m_NoiseSampRate = RATE_320;
        m_lpfir.Init( 320.0, m_Spread );
        m_LPGain = sqrt(320.0/(4.0*m_Spread*KGNB) );
    }
    else if( (m_Spread > 0.4) && (m_Spread <= 2.0) )
    {
        m_NoiseSampRate = RATE_64;
        m_lpfir.Init( 64.0, m_Spread );
        m_LPGain = sqrt(64.0/(4.0*m_Spread*KGNB) );
    }
    else if( (m_Spread >= 0.1) && (m_Spread <= 0.4) )
    {
        m_NoiseSampRate = RATE_12_8;
        m_lpfir.Init( 12.8, m_Spread );
        m_LPGain = sqrt(12.8/(4.0*m_Spread*KGNB) );
    }
    else if( (m_Spread >= 0.0) && (m_Spread < 0.1) )
    {		//here if spread<.1 so will not use any spread just offset
        m_NoiseSampRate = RATE_320;
        m_LPGain = 1.0;
    }
    for(int i=0; i<INTP_QUE_SIZE; i++)
    {
        m_pQue0[i].x = 0.0; m_pQue0[i].y = 0.0;
        m_pQue1[i].x = 0.0; m_pQue1[i].y = 0.0;
        m_pQue2[i].x = 0.0; m_pQue2[i].y = 0.0;
        m_pQue3[i].x = 0.0; m_pQue3[i].y = 0.0;
    }
    m_LPGain = m_LPGain/ sqrt((double)numpaths);
    for(int i=0; i<250; i++)
        MakeGaussianDelaySample();		//pre load filter
}

//////////////////////////////////////////////////////////////////////
// Performs a path calculation on pIn and puts it in pOut
//
//  Two Low Pass filtered Gaussian random numbers are created at
//	12.8, 64 Hz, or 320 Hz rate.  These form the input to a complex
//	interpolation filter that bumps the sample rate up to 8000Hz.
//
//	Two, three, or four stages of X5 upsampling/interpolation are used.
//	The complex noise is then multiplied by the input I/Q signal
//	to produce the spreading/fading simulation.
//
//  Finally a complex NCO is multiplied by the signal to produce a
//	Frequency offset.
//////////////////////////////////////////////////////////////////////
void Path::CalcPath(cmplx *pIn, cmplx *pOut)
{
int i,j;
cmplx acc;
cmplx tmp;
const double* Kptr;
cmplx* Firptr;
cmplx offset;
    if(m_PathActive)		// if this path is active
    {
        for(i=0; i<m_BlockSize; i++)
        {
            if( m_NoiseSampRate == RATE_12_8)
            {
                if( m_Indx%(5*5*5*5) == 0 )
                {			//generate noise samples at 12.8Hz rate
                    acc = MakeGaussianDelaySample();

//SweepGenCpx(  &acc, 12.8, 0.0, 6.4, 0.016 );

                    j = m_FirState0/INTP_VALUE;
                    m_pQue0[j].x = acc.x;
                    m_pQue0[j].y = acc.y;
                }
            }
            if( m_NoiseSampRate <= RATE_64)
            {
                if( m_Indx%(5*5*5) == 0 )
                {
                    if( m_NoiseSampRate == RATE_64)
                    {			//generate noise samples at 64Hz rate
                        acc = MakeGaussianDelaySample();
                    }
                    else
                    {
                        acc.x = 0.0; acc.y = 0.0;
                        Firptr = m_pQue0;
                        Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState0;
                        for(j=0; j<INTP_QUE_SIZE; j++)
                        {
                            acc.x += ( (Firptr->x)*(*Kptr) );
                            acc.y += ( (Firptr++->y)*(*Kptr) );
                            Kptr += INTP_VALUE;
                        }
                        if( --m_FirState0 < 0)
                            m_FirState0 = INTP_FIR_SIZE-1;
                    }

//SweepGenCpx(  &acc, 64, 0.0, 32.0, 0.08 );

                    j = m_FirState1/INTP_VALUE;
                    m_pQue1[j].x = acc.x;
                    m_pQue1[j].y = acc.y;
                }
            }
            if( m_Indx%(5*5) == 0 )	//interpolate/upsample x5
            {
                if( m_NoiseSampRate == RATE_320)
                {
                    acc = MakeGaussianDelaySample();
                }
                else
                {
                        acc.x = 0.0; acc.y = 0.0;
                        Firptr = m_pQue1;
                        Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState1;
                        for(j=0; j<INTP_QUE_SIZE; j++)
                        {
                            acc.x += ( (Firptr->x)*(*Kptr) );
                            acc.y += ( (Firptr++->y)*(*Kptr) );
                            Kptr += INTP_VALUE;
                        }
                        if( --m_FirState1 < 0)
                            m_FirState1 = INTP_FIR_SIZE-1;
                }

//SweepGenCpx(  &acc, 320, 0.0, 160.0, 0.4 );

                j = m_FirState2/INTP_VALUE;
                m_pQue2[j].x = acc.x;
                m_pQue2[j].y = acc.y;
            }
            if( m_Indx%(5) == 0 )	//interpolate/upsample x5
            {
                acc.x = 0.0; acc.y = 0.0;
                Firptr = m_pQue2;
                Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState2;
                for(j=0; j<INTP_QUE_SIZE; j++)
                {
                    acc.x += ( (Firptr->x)*(*Kptr) );
                    acc.y += ( (Firptr++->y)*(*Kptr) );
                    Kptr += INTP_VALUE;
                }
                if( --m_FirState2 < 0)
                    m_FirState2 = INTP_FIR_SIZE-1;

//SweepGenCpx(  &acc, 1600, 0.0, 800.0, 2 );

                j = m_FirState3/INTP_VALUE;
                m_pQue3[j].x = acc.x;
                m_pQue3[j].y = acc.y;
            }
            acc.x = 0.0; acc.y = 0.0;
            Firptr = m_pQue3;
            Kptr = X5IntrpFIRCoef+INTP_FIR_SIZE-m_FirState3;
            for(j=0; j<INTP_QUE_SIZE; j++)
            {
                acc.x += ( (Firptr->x)*(*Kptr) );
                acc.y += ( (Firptr++->y)*(*Kptr) );
                Kptr += INTP_VALUE;
            }
            if( --m_FirState3 < 0)
                m_FirState3 = INTP_FIR_SIZE-1;

//CalcCpxSweepRMS( acc, 8000);

            tmp.x = (acc.x*pIn[i].x - acc.y*pIn[i].y);
            tmp.y = (acc.x*pIn[i].y + acc.y*pIn[i].x);
            offset.x = cos(m_Timeinc);		//Cpx multiply by offset frequency
            offset.y = sin(m_Timeinc);
            pOut[i].x = ((offset.x*tmp.x) - (offset.y*tmp.y));
            pOut[i].y = ((offset.x*tmp.y) + (offset.y*tmp.x));
            m_Timeinc += (OFFSET_FREQ_CONST*m_Offset);
            m_Timeinc = fmod(m_Timeinc,K_2PI);	//keep radian counter bounded
            if( ++m_Indx > (INTP_VALUE*INTP_VALUE*INTP_VALUE*INTP_VALUE*m_BlockSize) )
                m_Indx = 0;
        }
    }
    else		// if path is not active just zero the output
    {
        for(i=0; i<m_BlockSize; i++)
        {
            pOut[i].x = 0.0;
            pOut[i].y = 0.0;
        }
    }
}

// Create the complex Rayleigh distributed samples by
// creating two Gaussian random distributed numbers for the I and Q
// terms and then passing them through a Gaussian shaped LP IIR.
// The 2 Sigma bandwidth of the LP filter determines the amount of spread.
cmplx Path::MakeGaussianDelaySample()
{
    cmplx val;
    if (m_Spread >= 0.1) {
        // Generate two uniform random numbers between -1 and +1 that are inside the unit circle
        double r2;
        do {
            val.x = 1.0 - 2.0 * (double)rand()/(double)RAND_MAX;
            val.y = 1.0 - 2.0 * (double)rand()/(double)RAND_MAX;
            r2 = val.x * val.x + val.y * val.y;
        } while (r2 >= 1.0 || r2 == 0.0);

        double scale = m_LPGain * sqrt(- 2.0 * log(r2) / r2);
        val.x *= scale;
        val.y *= scale;

        //SweepGenCpx(  &val, 320, 0.0, 30*5, 30*5/200.0);

        // Now LP filter the Gaussian samples
        val = m_lpfir.CalcFilter(val);
    } else
    {
        // Not using any spread.
        val.x = m_LPGain;
        val.y = 0;
    }

    //gDebug1 = CalcCpxRMS( val, 288000);
    //CalcCpxSweepRMS( val, 500);
    return val;
}

} // namespace PathSim
