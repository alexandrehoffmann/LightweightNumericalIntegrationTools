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
        const Scalar I   = Scalar( std::reduce(m_subIntergrals.cbegin(),    m_subIntergrals.cend(),    LongScalar(0)) );
        const Scalar err = Scalar( std::reduce(m_subIntergralsErr.cbegin(), m_subIntergralsErr.cend(), LongScalar(0)) );

        if (m_out) { fmt::print(m_out, "{} {:10.4e} {:10.4e} {:10.4e}\n", m_it, I, err, std::abs(I)*m_tol); }

        if (not std::isfinite(I))    { return I; }
        if (err < std::abs(I)*m_tol) { m_hasConverged = true; return I; }

        // we find the interval over which the integral is the least accurate
		const const_Iterator maxErrIt = std::max_element(m_subIntergralsErr.cbegin(), m_subIntergralsErr.cend());
		const Size maxErrIdx = Size(std::distance(m_subIntergralsErr.cbegin(), maxErrIt));
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
    return Scalar(std::reduce(m_subIntergrals.begin(), m_subIntergrals.end(), LongScalar(0)));
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateLeftInfinite(const Function& f, const Scalar xmax) -> Scalar
{
	using NumTraits = std::numeric_limits<Scalar>;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	if (std::abs(f(-NumTraits::infinity())) > NumTraits::epsilon())              { return NumTraits::quiet_NaN(); }	
	if (std::abs(gLaguerreQuad.integrateLeftInfinite(f)) < NumTraits::epsilon()) { return integrate(f, 0, xmax);  }
	
	Scalar xmin = -1;
	while (std::abs(gLaguerreQuad.integrateLeftInfinite(f, xmin)) >= NumTraits::epsilon())
	{
		xmin *= 2;
	}
	return integrate(f, xmin, xmax);
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrateRightInfinite(const Function& f, const Scalar xmin) -> Scalar
{
	using NumTraits = std::numeric_limits<Scalar>;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	if (std::abs(f(NumTraits::infinity())) > NumTraits::epsilon())                { return NumTraits::quiet_NaN(); }	
	if (std::abs(gLaguerreQuad.integrateRightInfinite(f)) < NumTraits::epsilon()) { return integrate(f, xmin, 0);  }
	
	Scalar xmax = 1;
	while (std::abs(gLaguerreQuad.integrateRightInfinite(f, xmax)) >= NumTraits::epsilon())
	{
		xmax *= 2;
	}
	return integrate(f, xmin, xmax);
}

template<class Derived> template<class Function>
auto AdaptiveQuadratureBase<Derived>::integrate(const Function& f) -> Scalar
{
	using NumTraits = std::numeric_limits<Scalar>;
	
	GaussLaguerreQuadrature<Scalar,LongScalar> gLaguerreQuad;
	
	if (std::abs(f(-NumTraits::infinity())) > NumTraits::epsilon())              { return NumTraits::quiet_NaN(); }	
	if (std::abs(f( NumTraits::infinity())) > NumTraits::epsilon())              { return NumTraits::quiet_NaN(); }	
	
	if (std::abs(gLaguerreQuad.integrateLeftInfinite(f)) < NumTraits::epsilon()) 
	{
		return std::abs(gLaguerreQuad.integrateRightInfinite(f)) < NumTraits::epsilon() ? 0 : integrateRightInfinite(f, 0); 
	}
	
	Scalar xmin = -1;
	while (std::abs(gLaguerreQuad.integrateLeftInfinite(f, xmin)) >= NumTraits::epsilon())
	{
		xmin *= 2;
	}
	
	return integrateRightInfinite(f, xmin);
}

} // namespace LNIT

#endif // LNIT_ADAPTIVE_QUADRATURE_BASE_IMPL_HPP
