#ifndef LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisAdaptiveQuadrature.hpp>

#include <LNIT/misc/StridedView.hpp>

namespace LNIT
{

//// explicit template instanciations ////
extern template class ClenshawCurtisAdaptiveQuadrature<float, float>;
extern template class ClenshawCurtisAdaptiveQuadrature<float, double>;
extern template class ClenshawCurtisAdaptiveQuadrature<float, long double>;

extern template class ClenshawCurtisAdaptiveQuadrature<double, double>;
extern template class ClenshawCurtisAdaptiveQuadrature<double, long double>;

extern template class ClenshawCurtisAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
auto ClenshawCurtisAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{	
    const auto fx = s_xi | std::views::transform([&f, xmin, xmax](const Scalar xi) -> Scalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	std::ranges::copy(fx, std::begin(m_fx33));
    
#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 202302L) || __cplusplus >= 202302L)
    const auto fx17 = m_fx33 | std::views::stride(2);
    const auto fx09 = m_fx33 | std::views::stride(4);
#else
    const misc::StridedView fx17(m_fx33, 2);
    const misc::StridedView fx09(m_fx33, 4);
#endif // if using c++23

    const LongScalar I33 = LongScalar(0.5)*(xmax - xmin)*std::inner_product(std::cbegin(s_wi33), std::cend(s_wi33), std::cbegin(m_fx33), LongScalar(0));
    const LongScalar I17 = LongScalar(0.5)*(xmax - xmin)*std::inner_product(std::cbegin(s_wi17), std::cend(s_wi17), std::cbegin(  fx17), LongScalar(0));
    const LongScalar I09 = LongScalar(0.5)*(xmax - xmin)*std::inner_product(std::cbegin(s_wi09), std::cend(s_wi09), std::cbegin(  fx09), LongScalar(0));
	
    const LongScalar err1 = std::abs(I33 - I17);
	const LongScalar err2 = std::abs(I33 - I09);
	
    return std::make_pair(I33, err2 == 0 ? LongScalar(0) : err1*(err1 / err2)*(err1 / err2));
    //~ return std::make_pair(I33, err2 == 0 ? LongScalar(0) : err2*(err1 / err2)*(err1 / err2));
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
