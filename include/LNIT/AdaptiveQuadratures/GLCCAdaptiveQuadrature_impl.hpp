#ifndef LNIT_GL_CC_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_GL_CC_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/GLCCAdaptiveQuadrature.hpp>

namespace LNIT
{

//// explicit template instanciations ////

extern template class GLCCAdaptiveQuadrature<float, float>;
extern template class GLCCAdaptiveQuadrature<float, double>;
extern template class GLCCAdaptiveQuadrature<float, long double>;

extern template class GLCCAdaptiveQuadrature<double, double>;
extern template class GLCCAdaptiveQuadrature<double, long double>;

extern template class GLCCAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
constexpr auto GLCCAdaptiveQuadrature<T,TT>::estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax) const -> std::pair<LongScalar, LongScalar>
{
    const auto fx_gl = s_xi_gl | std::views::transform([&f, xmin, xmax](const Scalar xi) -> Scalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	
    const auto fx_cc = s_xi_cc | std::views::transform([&f, xmin, xmax](const Scalar xi) -> Scalar
	{
		const Scalar x = Scalar(0.5)*(xi*(xmax - xmin) + (xmax + xmin));
		return f(x); 
	});
	
	const LongScalar I_gl = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_wi_gl), std::ranges::end(s_wi_gl), std::ranges::begin(fx_gl), LongScalar{});
	const LongScalar I_cc = LongScalar(0.5)*LongScalar(xmax - xmin)*std::inner_product(std::ranges::begin(s_wi_cc), std::ranges::end(s_wi_cc), std::ranges::begin(fx_cc), LongScalar{});
	
	return std::make_pair(I_gl, std::abs(I_cc - I_gl));
}

} // namespace LNIT

#endif // LNIT_GL_CC_ADAPTIVE_QUADRATURE_IMPL_HPP

