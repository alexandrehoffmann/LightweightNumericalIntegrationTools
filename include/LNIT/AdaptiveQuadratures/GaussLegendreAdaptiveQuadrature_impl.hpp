#ifndef LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GaussLegendreAdaptiveQuadrature<float, float>;
extern template class GaussLegendreAdaptiveQuadrature<float, double>;
extern template class GaussLegendreAdaptiveQuadrature<float, long double>;

extern template class GaussLegendreAdaptiveQuadrature<double, double>;
extern template class GaussLegendreAdaptiveQuadrature<double, long double>;

extern template class GaussLegendreAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
auto GaussLegendreAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{
    const auto fx = s_xi | std::views::transform([&f, xmin, xmax](const Scalar xi) -> Scalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	std::ranges::copy(fx, std::begin(m_fx15));
	
	const LongScalar I15 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::cbegin(s_wi15), std::cend(s_wi15), std::cbegin(m_fx15), LongScalar());
	const LongScalar I14 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::cbegin(s_wi14), std::cend(s_wi14), std::cbegin(m_fx15), LongScalar());
	const LongScalar I06 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::cbegin(s_wi06), std::cend(s_wi06), std::cbegin(m_fx15), LongScalar());

    const LongScalar err1 = std::abs(I15 - I14);
	const LongScalar err2 = std::abs(I15 - I06);

    return std::make_pair(I15, err2 == 0 ? LongScalar() : err1*(err1 / err2)*(err1 / err2));
    //~ return std::make_pair(I15, err2 == 0 ? LongScalar() : err2*(err1 / err2)*(err1 / err2));
}

} // namespace LNIT

#endif // LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_IMPL_HPP

