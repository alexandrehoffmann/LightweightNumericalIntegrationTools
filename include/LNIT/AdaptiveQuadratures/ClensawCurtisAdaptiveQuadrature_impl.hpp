#ifndef LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/ClensawCurtisAdaptiveQuadrature.hpp>

#include <LNIT/misc/StridedView.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class ClensawCurtisAdaptiveQuadrature<double, double>;
extern template class ClensawCurtisAdaptiveQuadrature<double, long double>;

extern template class ClensawCurtisAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
auto ClensawCurtisAdaptiveQuadrature<T,TT>::estimateIntegral_impl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{	
    for (Size i=0; i!=s_xi.size(); ++i)
    {
        const Scalar x = 0.5*(s_xi[i]*(xmax - xmin) + (xmax + xmin));
        m_fx33[i] = f(x);
    }
    
    const misc::StridedView fx17(m_fx33, 2);
    const misc::StridedView fx09(m_fx33, 4);
	
    const LongScalar I33 = 0.5*(xmax - xmin)*std::inner_product(s_wi33.cbegin(), s_wi33.cend(), m_fx33.cbegin(), LongScalar(0));
    const LongScalar I17 = 0.5*(xmax - xmin)*std::inner_product(s_wi17.cbegin(), s_wi17.cend(),   fx17.cbegin(), LongScalar(0));
    const LongScalar I09 = 0.5*(xmax - xmin)*std::inner_product(s_wi09.cbegin(), s_wi09.cend(),   fx09.cbegin(), LongScalar(0));
	
    const LongScalar err1 = std::abs(I33 - I17);
	const LongScalar err2 = std::abs(I33 - I09);
	
    return std::make_pair(I33, err2 == 0 ? LongScalar(0) : err2*(err1 / err2));
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
