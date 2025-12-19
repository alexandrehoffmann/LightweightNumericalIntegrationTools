#ifndef LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP
#define LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP

#include <cmath>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{

template<typename T, typename TT> class ClenshawCurtisAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< ClenshawCurtisAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

/**
 * @brief Adaptive Clenshaw–Curtis quadrature for numerical integration.
 *
 * This class implements adaptive quadrature using Clenshaw–Curtis rules
 * with nested evaluation points (9, 17, and 33 point formulas).
 *
 * @tparam T Floating point type for integration (e.g., double).
 * @tparam TT Higher precision type for accumulation (e.g., long double).
 */
template<typename T, typename TT=T>
class ClenshawCurtisAdaptiveQuadrature : public AdaptiveQuadratureBase< ClenshawCurtisAdaptiveQuadrature<T,TT> >
{
	using Base = AdaptiveQuadratureBase< ClenshawCurtisAdaptiveQuadrature<T,TT> >;
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
	template<class Function> std::pair<LongScalar, LongScalar> estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax);
	
	Scalar getMaxDeltaXImpl(const Scalar xmin, const Scalar xmax) const { return (xmax - xmin)*misc::maxDiff(std::span{s_xi}); } 
private:
	std::array<Scalar, 33> m_fx33;

	static constexpr std::array<Scalar, 9>  s_wi09 = {
		0.013888888888888888889, 0.15018690318427212327, 0.27539682539682539683, 0.36568611268874374975, 0.38968253968253968254, 
		0.36568611268874374975, 0.27539682539682539683, 0.15018690318427212327, 0.013888888888888888889};
	
	static constexpr std::array<Scalar, 17> s_wi17 = {
		0.0036764705882352941176, 0.037858898915636982870, 0.074992135464720461864, 0.10939574866032230299, 0.13846627228980170157, 
		0.16366286036013467511, 0.18098358815806198445, 0.19300406069135701942, 0.19591992974345915522, 0.19300406069135701942, 
		0.18098358815806198445, 0.16366286036013467511, 0.13846627228980170157, 0.10939574866032230299, 0.074992135464720461864, 
		0.037858898915636982870, 0.0036764705882352941176};
	
	static constexpr std::array<Scalar, 33> s_wi33 = {
		0.00094696969696969696970, 0.0094542927821143499895, 0.019173150313521813895, 0.028519011496393025382, 0.037533247094887870728, 
		0.046323857656911084780, 0.054493921484820975149, 0.062333204364453339841, 0.069366480811276115611, 0.075944895260547805768, 
		0.081573722835779174940, 0.086638633260986770727, 0.090645018048561663449, 0.094004339257928071018, 0.096231231126328843905, 
		0.097759283027214917108, 0.098117482962608961479, 0.097759283027214917108, 0.096231231126328843905, 0.094004339257928071018, 
		0.090645018048561663449, 0.086638633260986770727, 0.081573722835779174940, 0.075944895260547805768, 0.069366480811276115611, 
		0.062333204364453339841, 0.054493921484820975149, 0.046323857656911084780, 0.037533247094887870728, 0.028519011496393025382, 
		0.019173150313521813895, 0.0094542927821143499895, 0.00094696969696969696970};
		
	static constexpr std::array<Scalar, 33> s_xi = {
		1.0000000000000000000, 0.99518472667219688624, 0.98078528040323044913, 0.95694033573220886494, 0.92387953251128675613, 
		0.88192126434835502971, 0.83146961230254523708, 0.77301045336273696081, 0.70710678118654752440, 0.63439328416364549822, 
		0.55557023301960222474, 0.47139673682599764856, 0.38268343236508977173, 0.29028467725446236764, 0.19509032201612826785, 
		0.098017140329560601994, 0, -0.098017140329560601994, -0.19509032201612826785, -0.29028467725446236764, -0.38268343236508977173, 
		-0.47139673682599764856, -0.55557023301960222474, -0.63439328416364549822, -0.70710678118654752440, -0.77301045336273696081, 
		-0.83146961230254523708, -0.88192126434835502971, -0.92387953251128675613, -0.95694033573220886494, -0.98078528040323044913, 
		-0.99518472667219688624, -1.0000000000000000000};
};
	
} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClenshawCurtisAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP
