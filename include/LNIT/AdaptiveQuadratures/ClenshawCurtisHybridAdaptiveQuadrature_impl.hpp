#ifndef LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisHybridAdaptiveQuadrature.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class ClenshawCurtisHybridAdaptiveQuadrature<float, float>;
extern template class ClenshawCurtisHybridAdaptiveQuadrature<float, double>;
extern template class ClenshawCurtisHybridAdaptiveQuadrature<float, long double>;

extern template class ClenshawCurtisHybridAdaptiveQuadrature<double, double>;
extern template class ClenshawCurtisHybridAdaptiveQuadrature<double, long double>;

extern template class ClenshawCurtisHybridAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
auto ClenshawCurtisHybridAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{	
    const auto fx = s_xi | std::views::transform([&f, xmin, xmax](const Scalar xi) -> Scalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	std::ranges::copy(fx, std::begin(m_fx));
    
    const LongScalar I1 = Scalar(0.5)*(xmax - xmin)*std::inner_product(std::cbegin(s_wi),          std::cend(s_wi),          std::cbegin(m_fx), LongScalar(0));
    const LongScalar I2 = Scalar(0.5)*(xmax - xmin)*std::inner_product(std::cbegin(s_alternateWi), std::cend(s_alternateWi), std::cbegin(m_fx), LongScalar(0));
	
    return std::make_pair(I1, std::abs(I1 - I2));
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP
