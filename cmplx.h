#ifndef PATHSIM_CMPLX_HPP
#define PATHSIM_CMPLX_HPP

namespace PathSim {

struct cmplx {
	double r { 0. };
	double i { 0. };

	void   set(double r, double i) { this->r = r; this->i = i; }
	double l2() const { return r * r + i * i; }

	cmplx& operator+=(const cmplx &rhs) {
		this->r += rhs.r;
		this->i += rhs.i;
		return *this;
	}

	cmplx& operator*=(const double rhs) {
		this->r *= rhs;
		this->i *= rhs;
		return *this;
	}
};

static inline cmplx operator*(const cmplx &lhs, const cmplx &rhs)
{
	return cmplx{
        lhs.r * rhs.r - lhs.i * rhs.i,
        lhs.r * rhs.i + lhs.i * rhs.r
    };
}

static inline cmplx operator*(const cmplx &lhs, const double rhs)
{
	return cmplx{ lhs.r * rhs, lhs.i * rhs };
}

} // namespace PathSim

#endif // PATHSIM_CMPLX_HPP
