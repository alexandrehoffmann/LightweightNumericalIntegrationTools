#ifndef LNIT_GAUSS_LAGUERRE_QUADRATURE_IMPL_HPP
#define LNIT_GAUSS_LAGUERRE_QUADRATURE_IMPL_HPP

#include <numeric>
#include <fmt/core.h>
#include <ranges>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GaussLaguerreQuadrature<double, double>;
extern template class GaussLaguerreQuadrature<double, long double>;

extern template class GaussLaguerreQuadrature<long double, long double>;

//// method implementations ////

template<typename Scalar, typename LongScalar> template<class Function> 
Scalar GaussLaguerreQuadrature<Scalar, LongScalar>::integrateLeftInfinite(const Function& f, const Scalar a) const
{
	const auto fx = s_xi | std::views::transform([&f, a](const Scalar x) -> Scalar
	{
		return f(a - x); 
	});

	return Scalar(std::inner_product(std::cbegin(s_wi), std::cend(s_wi), std::cbegin(fx), LongScalar{}));	
}

template<typename Scalar, typename LongScalar> template<class Function> 
Scalar GaussLaguerreQuadrature<Scalar, LongScalar>::integrateRightInfinite(const Function& f, const Scalar a) const
{
	const auto fx = s_xi | std::views::transform([&f, a](const Scalar x) -> Scalar
	{
		return f(x + a); 
	});

	return Scalar(std::inner_product(std::cbegin(s_wi), std::cend(s_wi), std::cbegin(fx), LongScalar{}));	
}

} // namespace LNIT

#endif // LNIT_GAUSS_LAGUERRE_QUADRATURE_IMPL_HPP
