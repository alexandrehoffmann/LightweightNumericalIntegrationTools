#ifndef LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP
#define LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP

#include <cmath>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>

namespace LNIT
{

template<typename T, typename TT> class ClensawCurtisAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< ClensawCurtisAdaptiveQuadrature<T, TT> >
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
class ClensawCurtisAdaptiveQuadrature : public AdaptiveQuadratureBase< ClensawCurtisAdaptiveQuadrature<T,TT> >
{
	using Base = AdaptiveQuadratureBase< ClensawCurtisAdaptiveQuadrature<T,TT> >;
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
	template<class Function> std::pair<LongScalar, LongScalar> estimateIntegral_impl(const Function& f, const Scalar xmin, const Scalar xmax);
private:
	std::array<Scalar, 33> m_fx;

	static const std::array<Scalar, 9>  s_wi09;
	static const std::array<Scalar, 17> s_wi17;
	static const std::array<Scalar, 33> s_wi33;
	static const std::array<Scalar, 33> s_xi;
};
	
} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClensawCurtisAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP
