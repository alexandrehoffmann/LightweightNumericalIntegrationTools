#ifndef LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/ClensawCurtisHybridAdaptiveQuadrature.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class ClensawCurtisHybridAdaptiveQuadrature<double, double>;
extern template class ClensawCurtisHybridAdaptiveQuadrature<double, long double>;

extern template class ClensawCurtisHybridAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
auto ClensawCurtisHybridAdaptiveQuadrature<T,TT>::estimateIntegral_impl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{	
    for (Size i=0; i!=s_xi.size(); ++i)
    {
        const Scalar x = 0.5*(s_xi[i]*(xmax - xmin) + (xmax + xmin));
        m_fx[i] = f(x);
    }
    
    const LongScalar I1 = 0.5*(xmax - xmin)*std::inner_product(s_wi.cbegin(),          s_wi.cend(),          m_fx.cbegin(), LongScalar(0));
    const LongScalar I2 = 0.5*(xmax - xmin)*std::inner_product(s_alternateWi.cbegin(), s_alternateWi.cend(), m_fx.cbegin(), LongScalar(0));
	
    return std::make_pair(I1, std::abs(I1 - I2));
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_IMPL_HPP
