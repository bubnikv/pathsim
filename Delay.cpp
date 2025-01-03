// Delay.cpp: implementation of the CDelay class.
//   ( also performs Hilbert Real to complex I/Q 3KHz filtering )

#include "Delay.h"
#include "FilterTables.h"

namespace PathSim {

void Hilbert::init()
{
	memset(m_hilbert_queue, 0, sizeof(m_hilbert_queue));
	m_hilbert_ptr = HILBPFIR_LENGTH - 1;
}

// Hilbert 3KHz BP filters.  Real input and complex I/Q output
//   This FIR bandwidth limits the real input as well as creates a
//   complex I/Q output signal for the rest of the processing chain.
void Hilbert::filter_block(const double* pIn, cmplx* pOut)
{
	for (int i = 0; i < BLOCKSIZE; ++ i) {
		m_hilbert_queue[m_hilbert_ptr].set(pIn[i], pIn[i]);	//place real values in circular Queue
		const cmplx* Firptr = m_hilbert_queue;
		const double* IKptr = IHilbertBPFirCoef+HILBPFIR_LENGTH-m_hilbert_ptr;
		const double* QKptr = QHilbertBPFirCoef+HILBPFIR_LENGTH-m_hilbert_ptr;
		cmplx acc{0., 0.};
		for (int j = 0; j < HILBPFIR_LENGTH; ++ j, ++ Firptr)
			acc += (*Firptr) * (*IKptr++);
		pOut[i] = acc;
		if (-- m_hilbert_ptr < 0)
			m_hilbert_ptr = HILBPFIR_LENGTH - 1;
	}
}

void Delay::init()
{
	memset(m_delay_line, 0, sizeof(m_delay_line));
	m_in_ptr = BUFSIZE - 1;
	m_out_ptrs.clear();
}

void Delay::add_delay(double time_ms)
{
	m_out_ptrs.emplace_back(int(BUFSIZE - int(8.0 * time_ms) - 1));
}

// Uses pointers to create variable delays
void Delay::delay_block(const std::vector<cmplx> &inbuf, std::vector<std::vector<cmplx>*> &out_buffers)
{
    for (int i = 0; i < BLOCKSIZE; ++ i) {
		// Copy new data from inbuf into delay buffer
		m_delay_line[m_in_ptr ++] = inbuf[i];
		if (m_in_ptr >= BUFSIZE)
			m_in_ptr = 0;
		// Delay to the output buffers.
		for (int j = 0; j < m_out_ptrs.size(); ++ j) {
			int& ptr = m_out_ptrs[j];
			(*out_buffers[j])[i] = m_delay_line[ptr ++];
			if (ptr >= BUFSIZE)
				ptr = 0;
		}
	}
}

} // namespace PathSim
