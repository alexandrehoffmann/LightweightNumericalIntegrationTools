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
constexpr auto ClenshawCurtisHybridAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{	
	using std::abs;
	
    const auto fx = s_xi | std::views::transform([&f, xmin, xmax](const Scalar xi) -> Scalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	std::ranges::copy(fx, std::begin(m_fx));
    
    const LongScalar I1 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_wi),          std::ranges::end(s_wi),          std::ranges::begin(m_fx), LongScalar{});
    const LongScalar I2 = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_alternateWi), std::ranges::end(s_alternateWi), std::ranges::begin(m_fx), LongScalar{});

    return std::make_pair(I1, abs(I1 - I2));
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP
