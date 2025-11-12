#ifndef LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP
#define LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>

#include <cassert>
#include <numeric>
#include <algorithm>
#include <ranges>

#include <fmt/core.h>

#include <LNIT/GaussLaguerreQuadrature.hpp>

namespace LNIT
{

template<class Derived>
AdaptiveQuadratureBase<Derived>::AdaptiveQuadratureBase(const Size maxIt, const Scalar tol) 
    : m_maxIt(maxIt)
    , m_it(0)
    , m_tol(tol) 
{ 
    m_intervals.reserve(maxIt); 
    m_subIntergrals.reserve(maxIt); 
    m_subIntergralsErr.reserve(maxIt);
}

template<class Derived> template<class Function> 
auto AdaptiveQuadratureBase<Derived>::integrate(const Function& f, const Scalar xmin, const Scalar xmax) -> Scalar
{
    using const_Iterator = typename std::vector<LongScalar>::const_iterator;
    
    m_hasConverged = false;
	m_intervals.clear();
	m_subIntergrals.clear();
	m_subIntergralsErr.clear();

    Scalar res;
	Scalar estimatedErr;
	
	if (m_out) { fmt::print(m_out, "#NumericalIntegrator addapting quadrature over [{}, {}]\n", xmin, xmax); }
	if (m_out) { fmt::print(m_out, "#Iteration integral estimated_error tol\n"); }
	std::tie(res, estimatedErr) = estimateIntegral(f, xmin, xmax);
	
	m_intervals.emplace_back(xmin, xmax);
	m_subIntergrals.push_back(res);
	m_subIntergralsErr.push_back(estimatedErr);

    for (m_it=0; m_it!=m_maxIt; ++m_it)
	{
        const Scalar I   = Scalar( std::reduce(std::cbegin(m_subIntergrals),    std::cend(m_subIntergrals),    LongScalar(0)) );
        const Scalar err = Scalar( std::reduce(std::cbegin(m_subIntergralsErr), std::cend(m_subIntergralsErr), LongScalar(0)) );

        if (m_out) { fmt::print(m_out, "{} {:10.4e} {:10.4e} {:10.4e}\n", m_it, I, err, std::abs(I)*m_tol); }

        if (not std::isfinite(I))    { return I; }
        if (err < std::abs(I)*m_tol) { m_hasConverged = true; return I; }

        // we find the interval over which the integral is the least accurate
		const const_Iterator maxErrIt = std::max_element(std::cbegin(m_subIntergralsErr), std::cend(m_subIntergralsErr));
		const Size maxErrIdx = Size(std::distance(std::cbegin(m_subIntergralsErr), maxErrIt));
        // we split it in two
		const auto [a, b] = m_intervals[maxErrIdx];
		const Scalar midPoint = (a + b) / 2.;
        // first interval
		m_intervals[maxErrIdx] = Interval(a, midPoint);
		std::tie(m_subIntergrals[maxErrIdx], m_subIntergralsErr[maxErrIdx]) = estimateIntegral(f, a, midPoint);
        // second interval
		std::tie(res, estimatedErr) = estimateIntegral(f, midPoint, b);
		m_intervals.emplace_back(midPoint, b);
		m_subIntergrals.push_back(res);
		m_subIntergralsErr.push_back(estimatedErr);
    }
    return Scalar(std::reduce(std::cbegin(m_subIntergrals), std::cend(m_subIntergrals), LongScalar(0)));
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateLeftInfinite(const Function& f, const Scalar xmax) -> Scalar
{
	using NumTraits = std::numeric_limits<Scalar>;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	Scalar leftIntegral = gLaguerreQuad.integrateLeftInfinite(f);
	
	if (not std::isfinite(leftIntegral))               { return NumTraits::quiet_NaN(); }	
	if (std::abs(leftIntegral) < NumTraits::epsilon()) { return integrate(f, 0, xmax);  }
	
	Scalar xmin = -1;
	leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	while (std::isfinite(leftIntegral) and std::abs(leftIntegral) >= NumTraits::epsilon())
	{
		xmin *= 2;
		leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	}
	return std::isfinite(leftIntegral) ? integrate(f, xmin, xmax) : NumTraits::quiet_NaN();
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateRightInfinite(const Function& f, const Scalar xmin) -> Scalar
{
	using NumTraits = std::numeric_limits<Scalar>;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	Scalar rightIntegral = gLaguerreQuad.integrateRightInfinite(f);
	
	if (not std::isfinite(rightIntegral))               { return NumTraits::quiet_NaN(); }	
	if (std::abs(rightIntegral) < NumTraits::epsilon()) { return integrate(f, xmin, 0);  }
	
	Scalar xmax = 1;
	rightIntegral = gLaguerreQuad.integrateRightInfinite(f, xmax);
	while (std::isfinite(rightIntegral) and std::abs(rightIntegral) >= NumTraits::epsilon())
	{
		xmax *= 2;
		rightIntegral = gLaguerreQuad.integrateRightInfinite(f, xmax);
	}
	return std::isfinite(rightIntegral) ? integrate(f, xmin, xmax) : NumTraits::quiet_NaN();
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrate(const Function& f) -> Scalar
{
	using NumTraits = std::numeric_limits<Scalar>;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	Scalar leftIntegral = gLaguerreQuad.integrateLeftInfinite(f);
	
	if (not std::isfinite(leftIntegral)) { return NumTraits::quiet_NaN(); }	
	
	if (std::abs(leftIntegral) < NumTraits::epsilon()) 
	{
		const Scalar rightIntegral = gLaguerreQuad.integrateRightInfinite(f); 
		
		if (not std::isfinite(rightIntegral)) { return NumTraits::quiet_NaN(); } 
		return std::abs(rightIntegral) < NumTraits::epsilon() ? 0 : integrateRightInfinite(f, 0); 
	}
	
	Scalar xmin = -1;
	leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	while (std::isfinite(leftIntegral) and std::abs(leftIntegral) >= NumTraits::epsilon())
	{
		xmin *= 2;
		leftIntegral = gLaguerreQuad.integrateLeftInfinite(f, xmin);
	}
	
	return std::isfinite(leftIntegral) ? integrateRightInfinite(f, xmin) : NumTraits::quiet_NaN();
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::remapAndIntegrate(const Function& f) -> Scalar
{	
	const auto fref = [&f](const Scalar t) -> Scalar
	{
		return (t != -1 and t != 1) 
			? f(t / (1 - t*t))*(1 + t*t) / ((1 - t*t)*(1 - t*t))
			: 0;	
	};
	
	return integrate(fref, -1, 1);
}

} // namespace LNIT

#endif // LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP
