#ifndef LNIT_GL_CC_ADAPTIVE_QUADRATURE_HPP
#define LNIT_GL_CC_ADAPTIVE_QUADRATURE_HPP

#include <array>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{

template<typename T, typename TT> class GLCCAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< GLCCAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

/**
 * @class GLCCAdaptiveQuadrature
 * @brief Adaptive quadrature using combined Gauss–Legendre (GL) and Clenshaw–Curtis (CC) rules.
 *
 * This class computes two high-order estimates of the same integral — a
 * 15-point Gauss–Legendre and a 33-point Clenshaw–Curtis rule — which both are
 * of order 31. Because the two rules share the same order, their
 * difference provides a very accurate estimate of the integration error;
 * that error estimate is used to drive recursive adaptive subdivision.
 *
 * @tparam T  Scalar floating-point type (e.g., double).
 * @tparam TT Extended precision accumulator type (e.g., long double).
 */
template<typename T, typename TT=T> 
class GLCCAdaptiveQuadrature : public AdaptiveQuadratureBase< GLCCAdaptiveQuadrature<T, TT> >
{
	using Base = AdaptiveQuadratureBase< GLCCAdaptiveQuadrature<T,TT> >;
public:
	using Size       = Base::Size;       ///<  @brief Type for iteration counters.
	using Scalar     = Base::Scalar;     ///<  @brief Floating point type for integration (e.g., double).
	using LongScalar = Base::LongScalar; ///<  @brief Higher precision type for accumulation (e.g., long double).	
	
	/**
	 * @brief Estimate integral and error on [xmin, xmax].
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound.
	 * @param xmax Upper bound.
	 * @return Pair (integral, estimated error).
	 */
	template<class Function> std::pair<LongScalar, LongScalar> estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax) const;
	
	Scalar getMaxDeltaXImpl(const Scalar xmin, const Scalar xmax) const { return (xmax - xmin)*misc::maxDiff(std::span{s_xi_cc}); } 
private:
	static constexpr std::array<Scalar, 15> s_xi_gl{-0.98799251802048542849, -0.93727339240070590431, -0.84820658341042721620, -0.72441773136017004742, -0.57097217260853884754, -0.39415134707756336990, -0.20119409399743452230, 0, 0.20119409399743452230, 0.39415134707756336990, 0.57097217260853884754, 0.72441773136017004742, 0.84820658341042721620, 0.93727339240070590431, 0.98799251802048542849};
	static constexpr std::array<Scalar, 15> s_wi_gl{0.030753241996117268355, 0.070366047488108124709, 0.10715922046717193501, 0.13957067792615431445, 0.16626920581699393355, 0.18616100001556221103, 0.19843148532711157646, 0.20257824192556127288, 0.19843148532711157646, 0.18616100001556221103, 0.16626920581699393355, 0.13957067792615431445, 0.10715922046717193501, 0.070366047488108124709, 0.030753241996117268355};
	
	static constexpr std::array<Scalar, 33> s_xi_cc = {
		1.0000000000000000000, 0.99518472667219688624, 0.98078528040323044913, 0.95694033573220886494, 0.92387953251128675613, 
		0.88192126434835502971, 0.83146961230254523708, 0.77301045336273696081, 0.70710678118654752440, 0.63439328416364549822, 
		0.55557023301960222474, 0.47139673682599764856, 0.38268343236508977173, 0.29028467725446236764, 0.19509032201612826785, 
		0.098017140329560601994, 0, -0.098017140329560601994, -0.19509032201612826785, -0.29028467725446236764, -0.38268343236508977173, 
		-0.47139673682599764856, -0.55557023301960222474, -0.63439328416364549822, -0.70710678118654752440, -0.77301045336273696081, 
		-0.83146961230254523708, -0.88192126434835502971, -0.92387953251128675613, -0.95694033573220886494, -0.98078528040323044913, 
		-0.99518472667219688624, -1.0000000000000000000};
	static constexpr std::array<Scalar, 33> s_wi_cc = {
		0.00094696969696969696970, 0.0094542927821143499895, 0.019173150313521813895, 0.028519011496393025382, 0.037533247094887870728, 
		0.046323857656911084780, 0.054493921484820975149, 0.062333204364453339841, 0.069366480811276115611, 0.075944895260547805768, 
		0.081573722835779174940, 0.086638633260986770727, 0.090645018048561663449, 0.094004339257928071018, 0.096231231126328843905, 
		0.097759283027214917108, 0.098117482962608961479, 0.097759283027214917108, 0.096231231126328843905, 0.094004339257928071018, 
		0.090645018048561663449, 0.086638633260986770727, 0.081573722835779174940, 0.075944895260547805768, 0.069366480811276115611, 
		0.062333204364453339841, 0.054493921484820975149, 0.046323857656911084780, 0.037533247094887870728, 0.028519011496393025382, 
		0.019173150313521813895, 0.0094542927821143499895, 0.00094696969696969696970};

};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/GLCCAdaptiveQuadrature_impl.hpp>

#endif // LNIT_GL_CC_ADAPTIVE_QUADRATURE_HPP
