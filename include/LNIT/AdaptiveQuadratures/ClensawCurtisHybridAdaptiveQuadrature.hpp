#ifndef LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP
#define LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP

#include <cmath>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>

namespace LNIT
{

template<typename T, typename TT> class ClensawCurtisHybridAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< ClensawCurtisHybridAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

template<typename T, typename TT=T>
class ClensawCurtisHybridAdaptiveQuadrature : public AdaptiveQuadratureBase< ClensawCurtisHybridAdaptiveQuadrature<T,TT> >
{
	using Base = AdaptiveQuadratureBase< ClensawCurtisHybridAdaptiveQuadrature<T,TT> >;
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
	std::array<Scalar, 9> m_fx;

	static constexpr std::array<Scalar, 9> s_wi          = {1.0/42.0, 83.0/315.0, 47.0/105.0, 167.0/315.0, 47.0/105.0, 83.0/315.0, 1.0/42.0};
	static constexpr std::array<Scalar, 9> s_alternateWi = {13.0/45.0, -4.0/15.0, 44.0/45.0, 0, 44.0/45.0, -4.0/15.0, 13.0/45.0};
	static constexpr std::array<Scalar, 9> s_xi          = {1, (1.0/2.0)*std::sqrt(3), 1.0/2.0, 0, -1.0/2.0, -1.0/2.0*std::sqrt(3), -1};
};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClensawCurtisHybridAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP
