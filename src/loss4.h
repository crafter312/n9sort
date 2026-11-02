#ifndef CLoss4_H
#define CLoss4_H

#include <string>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>

class CLoss4 {
private:
	gsl_bspline_workspace* m_bw{nullptr};
	gsl_vector* m_coeffs{nullptr};

	size_t m_n{0}; // Number of basis splines (control points/coefficients)
	size_t m_k{0}; // Order (degree + 1)

public:
	CLoss4(const std::string&);
	~CLoss4();

	// Energy loss/gain functions
	double getEout(double, double, double);
	double getEin(double, double, double);

	// Disable copy operations to prevent double-free errors on GSL raw pointers
	CLoss4(const CLoss4&) = delete;
	CLoss4& operator=(const CLoss4&) = delete;

	// Allow move operations
	CLoss4(CLoss4&& other) noexcept 
		: m_bw(other.m_bw), m_coeffs(other.m_coeffs), m_n(other.m_n), m_k(other.m_k) {
		other.m_bw = nullptr;
		other.m_coeffs = nullptr;
	}
	CLoss4& operator=(CLoss4&& other) noexcept {
		if (this != &other) {
			if (m_bw) gsl_bspline_free(m_bw);
			if (m_coeffs) gsl_vector_free(m_coeffs);

			m_bw = other.m_bw;
			m_coeffs = other.m_coeffs;
			m_n = other.m_n;
			m_k = other.m_k;

			other.m_bw = nullptr;
			other.m_coeffs = nullptr;
		}
		return *this;
	}
};

#endif // CLoss4_H
