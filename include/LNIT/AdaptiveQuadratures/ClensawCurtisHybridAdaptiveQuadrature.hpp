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
	std::array<Scalar, 13> m_fx;

	static constexpr std::array<Scalar, 13> s_wi = {
		0.0064102564102564102564, 0.067222926117575560019, 0.13037703037703037703, 0.18592888592888592889, 0.22580752580752580753, 
		0.25384119494654550410, 0.26082436082436082436, 0.25384119494654550410, 0.22580752580752580753, 0.18592888592888592889, 
		0.13037703037703037703, 0.067222926117575560019, 0.0064102564102564102564};
	
	static constexpr std::array<Scalar, 13> s_alternateWi = {
		0, 0.096656546636198874454, 0.055555555555555555556, 0.32275132275132275132, 0.026984126984126984127, 0.49805244807279583454, 
		0, 0.49805244807279583454, 0.026984126984126984127, 0.32275132275132275132, 0.055555555555555555556, 0.096656546636198874454, 
		0};
		
	static constexpr std::array<Scalar, 13> s_xi = {
		1.0000000000000000000, 0.96592582628906828675, 0.86602540378443864676, 0.70710678118654752440, 0.50000000000000000000, 
		0.25881904510252076235, 0, -0.25881904510252076235, -0.50000000000000000000, -0.70710678118654752440, -0.86602540378443864676, 
		-0.96592582628906828675, -1.0000000000000000000};
};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClensawCurtisHybridAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_HYBRID_ADAPTIVE_QUADRATURE_HPP
