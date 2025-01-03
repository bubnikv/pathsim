#ifndef PATHSIM_DELAY_HPP
#define PATHSIM_DELAY_HPP

#include <math.h>
#include <vector>

#include "cmplx.h"
#include "FilterTables.h"

namespace PathSim {

class Hilbert
{
public:
    static constexpr int BLOCKSIZE  = 2048;

    void init();
    void filter_block(const double* pIn, cmplx* pOut);

private:
    cmplx m_hilbert_queue[HILBPFIR_LENGTH];
    int   m_hilbert_ptr = 0;
};

class Delay
{
public:
    // 50mSecs max delay
    static constexpr int MAXDELAY   = 50 * 8;
    static constexpr int BLOCKSIZE  = Hilbert::BLOCKSIZE;
    static constexpr int BUFSIZE    = BLOCKSIZE + MAXDELAY;

    void init();
    void add_delay(double time_ms);
    void delay_block(const std::vector<cmplx> &inbuf, std::vector<std::vector<cmplx>*> &out_buffers);

private:
    cmplx 	m_delay_line[BUFSIZE];
    int 	m_in_ptr = 0;
    std::vector<int> m_out_ptrs;
};

} // namespace PathSim

#endif // PATHSIM_DELAY_HPP
