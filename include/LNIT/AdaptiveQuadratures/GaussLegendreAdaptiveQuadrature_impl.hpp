#ifndef LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature.hpp>

#include <LNIT/misc/ArrayView.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GaussLegendreAdaptiveQuadrature<double, double>;
extern template class GaussLegendreAdaptiveQuadrature<double, long double>;

extern template class GaussLegendreAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
auto GaussLegendreAdaptiveQuadrature<T,TT>::estimateIntegral_impl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{
    const auto fx = s_xi | std::views::transform([&f, xmin, xmax](const Scalar xi) -> Scalar
	{
		const Scalar x = xmin + (xmax - xmin)*xi;
		return f(x); 
	});
	std::ranges::copy(fx, std::begin(m_fx15));
	
	const misc::ArrayView fx14(m_fx15, s_15To14);
	const misc::ArrayView fx06(m_fx15, s_15To06);
	
	const LongScalar I15 = (xmax - xmin)*std::inner_product(std::cbegin(s_wi15), std::cend(s_wi15), std::cbegin(m_fx15), LongScalar(0));
	const LongScalar I14 = (xmax - xmin)*std::inner_product(std::cbegin(s_wi14), std::cend(s_wi14), std::cbegin(  fx14), LongScalar(0));
	const LongScalar I06 = (xmax - xmin)*std::inner_product(std::cbegin(s_wi06), std::cend(s_wi06), std::cbegin(  fx06), LongScalar(0));

    const LongScalar err1 = std::abs(I15 - I14);
	const LongScalar err2 = std::abs(I15 - I06);

    return std::make_pair(I15, err2 == 0 ? LongScalar(0) : err1*(err1 / err2)*(err1 / err2));
    //~ return std::make_pair(I15, err2 == 0 ? LongScalar(0) : err2*(err1 / err2)*(err1 / err2));
}

} // namespace LNIT

#endif // LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP

