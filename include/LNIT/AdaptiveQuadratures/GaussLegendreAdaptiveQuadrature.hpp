#ifndef LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_HPP
#define LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_HPP

#include <array>

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase.hpp>
#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{

template<typename T, typename TT> class GaussLegendreAdaptiveQuadrature;

template<typename T, typename TT> 
struct AdaptiveQuadratureTraits< GaussLegendreAdaptiveQuadrature<T, TT> >
{
	using Size       = unsigned int;
	using Scalar     = T;
	using LongScalar = TT;
};

template<typename T, typename TT=T> 
class GaussLegendreAdaptiveQuadrature : public AdaptiveQuadratureBase< GaussLegendreAdaptiveQuadrature<T, TT> >
{
	using Base = AdaptiveQuadratureBase< GaussLegendreAdaptiveQuadrature<T,TT> >;
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
	std::array<LongScalar, 15> m_fx15;
	
	static constexpr std::array<Scalar, 15> s_xi{0.00600374099, 0.0313633038 , 0.075896708295, 0.13779113432, 0.214513913696, 0.302924326461, 0.399402953001, 0.5, 0.600597046999, 0.697075673539, 0.785486086304, 0.86220886568, 0.924103291705, 0.9686366962, 0.99399625901};
	
	static constexpr std::array<Scalar, 15> s_wi15{0.015376620998, 0.035183023744, 0.053579610234, 0.069785338963, 0.083134602908, 0.093080500008, 0.099215742664, 0.101289120963, 0.099215742664, 0.093080500008, 0.083134602908, 0.069785338963, 0.053579610234, 0.035183023744, 0.015376620998};
	static constexpr std::array<Scalar, 14> s_wi14{0.021474027016, 0.014373158877, 0.092599212015, 0.011827748921, 0.158470029591, 0.003842924025, 0.197412898929, 0.197412903702, 0.003842915366, 0.15847004048 , 0.011827737944, 0.092599220945, 0.014373153478, 0.02147402871};
	static constexpr std::array<Scalar, 6>  s_wi06{0.107157609486, 0.03113090193 , 0.361711488584, 0.361711488584, 0.03113090193 , 0.107157609486};	
	
	static constexpr std::array<Size, 14> s_15To14{0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14};
	static constexpr std::array<Size, 14> s_15To06{1, 3, 5, 9, 11, 13};
};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature_impl.hpp>

#endif // LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_HPP
