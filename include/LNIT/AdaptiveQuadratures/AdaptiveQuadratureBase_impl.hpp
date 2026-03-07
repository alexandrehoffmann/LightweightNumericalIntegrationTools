#ifndef LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP
#define LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/GaussLaguerreQuadrature.hpp>
#include <LNIT/misc/Numeric.hpp>

#include <cassert>
#include <algorithm>
#include <ranges>
#include <cmath>

#include <fmt/core.h>

namespace LNIT
{

template<class Derived>
AdaptiveQuadratureBase<Derived>::AdaptiveQuadratureBase(const Size maxIt, const Scalar relativeTol, const Scalar absoluteTol) 
    : m_maxIt(maxIt)
    , m_it(0)
    , m_relativeTol(relativeTol)
    , m_absoluteTol(absoluteTol) 
{ 
    m_intervals.reserve(maxIt); 
    m_subIntergrals.reserve(maxIt); 
    m_subIntergralsErr.reserve(maxIt);
}

template<class Derived> template<class Function> 
auto AdaptiveQuadratureBase<Derived>::integrate(const Function& f, const Scalar xmin, const Scalar xmax) -> LongScalar
{	
	using std::ceil;
	using std::abs;
	using std::isfinite;
	
    using const_Iterator = typename std::vector<LongScalar>::const_iterator;
    
	m_hasConverged = false;
	m_intervals.clear();
	m_subIntergrals.clear();
	m_subIntergralsErr.clear();

	LongScalar res;
	LongScalar estimatedErr;
	
	if (m_out) { fmt::print(m_out, "#NumericalIntegrator addapting quadrature over [{}, {}]\n", xmin, xmax); }
	if (m_out) { fmt::print(m_out, "#Iteration integral estimated_error relative_tol absolute_tol\n"); }
	
	const Size N = Size(ceil(getMaxDeltaX(xmin, xmax)));
	
	m_intervals.reserve(N);
	m_subIntergrals.reserve(N);
	m_subIntergralsErr.reserve(N);
	
	for (Size i=0; i!=N; ++i)
	{
		const Scalar x_i   = xmin + Scalar(i)*(xmax - xmin) / Scalar(N);
		const Scalar x_ip1 = xmin + Scalar(i+1)*(xmax - xmin) / Scalar(N);
		
		std::tie(res, estimatedErr) = estimateIntegral(f, x_i, x_ip1);
		
		m_intervals.emplace_back(x_i, x_ip1);
		m_subIntergrals.push_back(res);
		m_subIntergralsErr.push_back(estimatedErr);
	}
	
	for (m_it=0; m_it!=m_maxIt; ++m_it)
	{
		const LongScalar I   = getEstimatedIntegral();
		const LongScalar err = getEstimatedError();

		if (m_out) { fmt::print(m_out, "{} {:10.4e} {:10.4e} {:10.4e} {:10.4e}\n", m_it, Scalar(I), Scalar(err), Scalar(abs(I)*m_relativeTol), Scalar(m_absoluteTol)); }

		if (not isfinite(I)) { return I; }
		if (err < abs(I)*m_relativeTol or err < m_absoluteTol) { m_hasConverged = true; return I; }

		// we find the interval over which the integral is the least accurate
		const const_Iterator maxErrIt = std::max_element(std::ranges::cbegin(m_subIntergralsErr), std::ranges::cend(m_subIntergralsErr));
		const Size maxErrIdx = Size(std::ranges::distance(std::ranges::cbegin(m_subIntergralsErr), maxErrIt));
		// we split it in two
		const auto [a, b] = m_intervals[maxErrIdx];
		const Scalar midPoint = Scalar(0.5)*(a + b);
		// first interval
		m_intervals[maxErrIdx] = Interval(a, midPoint);
		std::tie(m_subIntergrals[maxErrIdx], m_subIntergralsErr[maxErrIdx]) = estimateIntegral(f, a, midPoint);
		// second interval
		std::tie(res, estimatedErr) = estimateIntegral(f, midPoint, b);
		m_intervals.emplace_back(midPoint, b);
		m_subIntergrals.push_back(res);
		m_subIntergralsErr.push_back(estimatedErr);
	}
	return getEstimatedIntegral();
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateLeftInfinite(const Function& f, const Scalar xmax) -> LongScalar
{
	using std::isfinite;
	using std::abs;
	
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	LongScalar leftIntegral = gLaguerreQuad.integrateLeftInfinite(f);
	
	if (not isfinite(leftIntegral))		                      { return NumTraits<LongScalar>::NaN; }	
	if (abs(leftIntegral) < NumTraits<LongScalar>::epsilon) { return integrate(f, Scalar{}, xmax);  }
	
	Scalar xmin = -1;
	leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	while (isfinite(leftIntegral) and abs(leftIntegral) >= NumTraits<LongScalar>::epsilon)
	{
		xmin *= 2;
		leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	}
	return isfinite(leftIntegral) ? integrate(f, xmin, xmax) : NumTraits<LongScalar>::NaN;
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateRightInfinite(const Function& f, const Scalar xmin) -> LongScalar
{
	using std::isfinite;
	using std::abs;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	LongScalar rightIntegral = gLaguerreQuad.integrateRightInfinite(f);
	
	if (not isfinite(rightIntegral))		                     { return NumTraits<LongScalar>::NaN; }	
	if (abs(rightIntegral) < NumTraits<LongScalar>::epsilon) { return integrate(f, xmin, Scalar{});  }
	
	Scalar xmax = 1;
	rightIntegral = gLaguerreQuad.integrateRightInfinite(f, xmax);
	while (isfinite(rightIntegral) and abs(rightIntegral) >= NumTraits<LongScalar>::epsilon)
	{
		xmax *= 2;
		rightIntegral = gLaguerreQuad.integrateRightInfinite(f, xmax);
	}
	return isfinite(rightIntegral) ? integrate(f, xmin, xmax) : NumTraits<LongScalar>::NaN;
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrate(const Function& f) -> LongScalar
{
	using std::isfinite;
	using std::abs;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	LongScalar leftIntegral = gLaguerreQuad.integrateLeftInfinite(f);
	
	if (not isfinite(leftIntegral)) { return NumTraits<LongScalar>::NaN; }	
	
	if (abs(leftIntegral) < NumTraits<LongScalar>::epsilon) 
	{
		const LongScalar rightIntegral = gLaguerreQuad.integrateRightInfinite(f); 
		
		if (not isfinite(rightIntegral)) { return NumTraits<LongScalar>::NaN; } 
		return abs(rightIntegral) < NumTraits<LongScalar>::epsilon ? LongScalar{} : integrateRightInfinite(f, Scalar{}); 
	}
	
	Scalar xmin = -1;
	leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	while (isfinite(leftIntegral) and abs(leftIntegral) >= NumTraits<LongScalar>::epsilon)
	{
		xmin *= 2;
		leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	}
	
	return isfinite(leftIntegral) ? integrateRightInfinite(f, xmin) : NumTraits<LongScalar>::NaN;
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::remapAndIntegrate(const Function& f) -> LongScalar
{		
	using std::isnan;
	
	const auto fref = [&f](const Scalar t) -> LongScalar
	{			
		const LongScalar fx = f(t / (1 - t*t));
			
		return isnan(fx)
			? LongScalar{}
			: fx*LongScalar(1 + t*t) / LongScalar((1 - t*t)*(1 - t*t));	
	};
	
	return integrate(fref, -1, 1);
}

} // namespace LNIT

#endif // LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP
