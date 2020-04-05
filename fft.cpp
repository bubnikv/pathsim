// fft.cpp: implementation of the Cfft class.
//  This is a slightly modified version of Takuya OOURA's
//     original radix 4 FFT package.
//Copyright(C) 1996-1998 Takuya OOURA
//    (email: ooura@mmm.t.u-tokyo.ac.jp).
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "fft.h"

//////////////////////////////////////////////////////////////////////
// Local Defines
//////////////////////////////////////////////////////////////////////
#define SQRT_FFT_SIZE 46//sqrt(2048)
#define K_2PI (8.0*atan(1))
#define FREQ_SCALE (8000.0/2048.0)
#define SAMPLE_RATE 8000

//////////////////////////////////////////////////////////////////////
// A pure input sin wave ... Asin(wt)... will produce an fft output 
//   peak of (N*A/4)^2  where N is FFT_SIZE.
// To convert to a Power dB range:
//   PdBmax = 10*log10( (N*A/4)^2 + K_C ) + K_B
//   PdBmin = 10*log10( 0 + K_C ) + K_B
//  if (N*A/4)^2 >> K_C 
//  Then K_B = PdBmax - 20*log10( N*A/4 )
//       K_C = 10 ^ ( (PdBmin-K_B)/10 )
//  for power range of 0 to 100 dB with input(A) of 32767 and N=2048
//			K_B = -44.494132  and K_C = 2.81458e4
// To eliminate the multiply by 10, divide by 10 so for an output
//		range of 0 to 100dB the stored value is 0.0 to 10.0
//   so final constant K_B = -4.4494132
#define K_B (-4.4494132)
#define K_C (2.81458e4)

#define K_ROOT (1.0/4.0)		
#define K_ROOTGN 289.626		//(A*N/8)^(2/4) / 10

namespace PathSim {

Cfft::Cfft()
{
	WindowTbl = new double[FFT_SIZE];
	SinCosTbl = new double[FFT_SIZE/2];
	WorkArea = new int[SQRT_FFT_SIZE+2];
	pFFTAveBuf = new double[FFT_SIZE];	//Even index's hold average
	pFFTInBuf = new double[FFT_SIZE];
	WorkArea[0] = 0;
	makewt(FFT_SIZE/4, WorkArea, SinCosTbl);
    makect(FFT_SIZE/4, WorkArea, SinCosTbl + WorkArea[0]);
	for(int i=0; i<FFT_SIZE; i++)
	{
		pFFTAveBuf[i] = 0.0;
// Pick a data windowing function:
//		WindowTbl[i] = 1.0;										//rectangle
//		WindowTbl[i] = .54 - .46*cos( (K_2PI*i)/(FFT_SIZE-1) );	//Hamming
		WindowTbl[i] = (.5 - .5*cos( (K_2PI*i)/(FFT_SIZE-1) ));	//Hanning
	}
	m_AveSize = 1;
	m_LastAve = 1;
	m_LogMode = true;
	m_LastLogMode = true;
	m_Gain = 10.0;
	m_Clip = 0.0;
	m_Overload = false;
}

Cfft::~Cfft()
{							// free all resources
	if(WorkArea)
	{
		delete WorkArea;
		WorkArea = NULL;
	}
	if(SinCosTbl)
	{
		delete SinCosTbl;
		SinCosTbl = NULL;
	}
	if(WindowTbl)
	{
		delete WindowTbl;
		WindowTbl = NULL;
	}
	if(pFFTAveBuf)
	{
		delete pFFTAveBuf;
		pFFTAveBuf = NULL;
	}
	if(pFFTInBuf)
	{
		delete pFFTInBuf;
		pFFTInBuf = NULL;
	}

}

void Cfft::SetFFTParams(int ave, double gain, int type )
{
	if( type==0 )
		m_LogMode = false;
	else
		m_LogMode = true;
	if( (type>=10) && (type<=99) )
		m_Clip = (double)type/10.0;
	else
		m_Clip = 0.0;
	if(ave>0)
		m_AveSize = ave;
	else
		m_AveSize = 1;
	if( (m_LastAve != ave) || (m_LastLogMode != m_LogMode) )
		ResetFFT();
	m_LastLogMode = m_LogMode;
	m_LastAve = m_AveSize;
	m_Gain = 0.1*(gain*10.0/(10.0-m_Clip));

}

void Cfft::ResetFFT()
{
	for(int i=0; i<FFT_SIZE;i++)
		pFFTAveBuf[i] = 0.0;

}


//////////////////////////////////////////////////////////////////////
// "InBuf[]" is first multiplied by a window function and then
//  calculates an "FFT_SIZE" point FFT on "InBuf[]".
//  The result is converted to dB or 4th root and stored in pFFTAveBuf
//  If "Ave" is > 1, a LP smoothing filter 
//  is calculated on the output.
//////////////////////////////////////////////////////////////////////
void Cfft::CalcFFT(double * InBuf)
{
int i;
	m_Overload = false;
	for(i=0; i<FFT_SIZE; i++)
	{
		if( InBuf[i] > 32768.0*0.90 )	//flag overload if within 10% of max
			m_Overload = true;
		pFFTInBuf[i] = WindowTbl[i] * InBuf[i];		//window the data
	}
//Calculate the FFT
	bitrv2(FFT_SIZE, WorkArea + 2, pFFTInBuf);
	cftfsub(FFT_SIZE, pFFTInBuf, SinCosTbl);
	rftfsub(FFT_SIZE, pFFTInBuf, WorkArea[1], SinCosTbl + WorkArea[0]);
}

//////////////////////////////////////////////////////////////////////
//	the range "start" to "stop" is multiplied by "gain" and copied
//   into "OutBuf[]".
// The function returns true if the input is overloaded
//////////////////////////////////////////////////////////////////////
bool Cfft::GetFFTData(int start, int stop, long* OutBuf )
{
	for( int i=start; i<=stop; i++ )		//copy and scale into OutBuf[]
		OutBuf[i] = (long)(m_Gain*pFFTAveBuf[i<<1]);
	return m_Overload;
}


// Nitty gritty fft routines by Takuya OOURA
void Cfft::rftfsub(int n, double *a, int nc, double *c)
{
double tmp;
    int j, k, kk, ks, m;
    double wkr, wki, xr, xi, yr, yi;
    
    ks = (nc << 2) / n;
    kk = 0;
    m = n >> 1;
	if(m_LogMode)
	{
		for (k = 2; k < m; k += 2 ) 
		{
			j = n - k;
			kk += ks;
			wkr = 0.5 - c[nc - kk];
			wki = c[kk];
			xr = a[k] - a[j];
			xi = a[k + 1] + a[j + 1];
			yr = wkr * xr - wki * xi;
			yi = wkr * xi + wki * xr;
			a[k] -= yr;
			xi = a[k]*a[k];
			a[k+1] -= yi;
			xi += ( a[k+1]*a[k+1]);
			a[j] += yr;
			xr = a[j]*a[j];
			a[j+1] -= yi;
			xr += (a[j+1]*a[j+1]);
			tmp = log10(xi+K_C) + K_B;
			if( (tmp -= m_Clip)<0.0 )
				pFFTAveBuf[k] = (1.0-1.0/m_AveSize)*pFFTAveBuf[k];
			else
 				pFFTAveBuf[k] = (1.0-1.0/m_AveSize)*pFFTAveBuf[k] +
									(1.0/m_AveSize)*tmp;
			tmp = log10(xr+K_C) + K_B;
			if( (tmp -= m_Clip)<0.0 )
				pFFTAveBuf[j] = (1.0-1.0/m_AveSize)*pFFTAveBuf[j];
			else
 				pFFTAveBuf[j] = (1.0-1.0/m_AveSize)*pFFTAveBuf[j] +
									(1.0/m_AveSize)*tmp;
		}
	}
	else
	{
		for (k = 2; k < m; k += 2 ) 
		{
			j = n - k;
			kk += ks;
			wkr = 0.5 - c[nc - kk];
			wki = c[kk];
			xr = a[k] - a[j];
			xi = a[k + 1] + a[j + 1];
			yr = wkr * xr - wki * xi;
			yi = wkr * xi + wki * xr;
			a[k] -= yr;
			xi = a[k]*a[k];
			a[k+1] -= yi;
			xi += ( a[k+1]*a[k+1]);
			a[j] += yr;
			xr = a[j]*a[j];
			a[j+1] -= yi;
			xr += (a[j+1]*a[j+1]);

			pFFTAveBuf[k] = (1.0-1.0/m_AveSize)*pFFTAveBuf[k] + 
					(1.0/m_AveSize)*pow(xi , K_ROOT)/K_ROOTGN;
 			pFFTAveBuf[j] = (1.0-1.0/m_AveSize)*pFFTAveBuf[j] +
					(1.0/m_AveSize)*pow(xr , K_ROOT)/K_ROOTGN;
		}
	}
	pFFTAveBuf[1024] = pFFTAveBuf[1022];
}

/* -------- initializing routines -------- */
void Cfft::makewt(int nw, int *ip, double *w)
{
    int nwh, j;
    double delta, x, y;
    
    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        for (j = 2; j < nwh; j += 2) {
            x = cos(delta * j);
            y = sin(delta * j);
            w[j] = x;
            w[j + 1] = y;
            w[nw - j] = y;
            w[nw - j + 1] = x;
        }
        bitrv2(nw, ip + 2, w);
    }
}


void Cfft::makect(int nc, int *ip, double *c)
{
    int nch, j;
    double delta;
    
    ip[1] = nc;
    if (nc > 1) {
        nch = nc >> 1;
        delta = atan(1.0) / nch;
        c[0] = cos(delta * nch);
        c[nch] = 0.5 * c[0];
        for (j = 1; j < nch; j++) {
            c[j] = 0.5 * cos(delta * j);
            c[nc - j] = 0.5 * sin(delta * j);
        }
    }
}


/* -------- child routines -------- */
void Cfft::bitrv2(int n, int *ip, double *a)
{
    int j, j1, k, k1, l, m, m2;
    double xr, xi;
    
    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 2) < l) {
        l >>= 1;
        for (j = 0; j < m; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    if ((m << 2) > l) {
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    } else {
        m2 = m << 1;
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = (j << 1) + ip[k];
                k1 = (k << 1) + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = a[j1 + 1];
                a[j1] = a[k1];
                a[j1 + 1] = a[k1 + 1];
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    }
}

void Cfft::cftfsub(int n, double *a, double *w)
{
    int j, j1, j2, j3, l;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i - x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
        }
    } else {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            x0r = a[j] - a[j1];
            x0i = a[j + 1] - a[j1 + 1];
            a[j] += a[j1];
            a[j + 1] += a[j1 + 1];
            a[j1] = x0r;
            a[j1 + 1] = x0i;
        }
    }
}



void Cfft::cft1st(int n, double *a, double *w)
{
    int j, k1, k2;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    x0r = a[0] + a[2];
    x0i = a[1] + a[3];
    x1r = a[0] - a[2];
    x1i = a[1] - a[3];
    x2r = a[4] + a[6];
    x2i = a[5] + a[7];
    x3r = a[4] - a[6];
    x3i = a[5] - a[7];
    a[0] = x0r + x2r;
    a[1] = x0i + x2i;
    a[4] = x0r - x2r;
    a[5] = x0i - x2i;
    a[2] = x1r - x3i;
    a[3] = x1i + x3r;
    a[6] = x1r + x3i;
    a[7] = x1i - x3r;
    wk1r = w[2];
    x0r = a[8] + a[10];
    x0i = a[9] + a[11];
    x1r = a[8] - a[10];
    x1i = a[9] - a[11];
    x2r = a[12] + a[14];
    x2i = a[13] + a[15];
    x3r = a[12] - a[14];
    x3i = a[13] - a[15];
    a[8] = x0r + x2r;
    a[9] = x0i + x2i;
    a[12] = x2i - x0i;
    a[13] = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[10] = wk1r * (x0r - x0i);
    a[11] = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[14] = wk1r * (x0i - x0r);
    a[15] = wk1r * (x0i + x0r);
    k1 = 0;
    for (j = 16; j < n; j += 16) {
        k1 += 2;
        k2 = k1 << 1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        x0r = a[j] + a[j + 2];
        x0i = a[j + 1] + a[j + 3];
        x1r = a[j] - a[j + 2];
        x1i = a[j + 1] - a[j + 3];
        x2r = a[j + 4] + a[j + 6];
        x2i = a[j + 5] + a[j + 7];
        x3r = a[j + 4] - a[j + 6];
        x3i = a[j + 5] - a[j + 7];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 4] = wk2r * x0r - wk2i * x0i;
        a[j + 5] = wk2r * x0i + wk2i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 2] = wk1r * x0r - wk1i * x0i;
        a[j + 3] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 6] = wk3r * x0r - wk3i * x0i;
        a[j + 7] = wk3r * x0i + wk3i * x0r;
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        x0r = a[j + 8] + a[j + 10];
        x0i = a[j + 9] + a[j + 11];
        x1r = a[j + 8] - a[j + 10];
        x1i = a[j + 9] - a[j + 11];
        x2r = a[j + 12] + a[j + 14];
        x2i = a[j + 13] + a[j + 15];
        x3r = a[j + 12] - a[j + 14];
        x3i = a[j + 13] - a[j + 15];
        a[j + 8] = x0r + x2r;
        a[j + 9] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 12] = -wk2i * x0r - wk2r * x0i;
        a[j + 13] = -wk2i * x0i + wk2r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 10] = wk1r * x0r - wk1i * x0i;
        a[j + 11] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 14] = wk3r * x0r - wk3i * x0i;
        a[j + 15] = wk3r * x0i + wk3i * x0r;
    }
}


void Cfft::cftmdl(int n, int l, double *a, double *w)
{
    int j, j1, j2, j3, k, k1, k2, m, m2;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    m = l << 2;
    for (j = 0; j < l; j += 2) {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x0r - x2r;
        a[j2 + 1] = x0i - x2i;
        a[j1] = x1r - x3i;
        a[j1 + 1] = x1i + x3r;
        a[j3] = x1r + x3i;
        a[j3 + 1] = x1i - x3r;
    }
    wk1r = w[2];
    for (j = m; j < l + m; j += 2) {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x2i - x0i;
        a[j2 + 1] = x0r - x2r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j1] = wk1r * (x0r - x0i);
        a[j1 + 1] = wk1r * (x0r + x0i);
        x0r = x3i + x1r;
        x0i = x3r - x1i;
        a[j3] = wk1r * (x0i - x0r);
        a[j3 + 1] = wk1r * (x0i + x0r);
    }
    k1 = 0;
    m2 = m << 1;
    for (k = m2; k < n; k += m2) {
        k1 += 2;
        k2 = k1 << 1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        for (j = k; j < l + k; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = wk2r * x0r - wk2i * x0i;
            a[j2 + 1] = wk2r * x0i + wk2i * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        for (j = k + m; j < l + (k + m); j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = -wk2i * x0r - wk2r * x0i;
            a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
    }
}

} // namespace PathSim
