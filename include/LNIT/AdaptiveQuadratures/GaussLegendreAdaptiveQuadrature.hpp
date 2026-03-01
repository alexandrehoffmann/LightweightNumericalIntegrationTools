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
	template<class Function> constexpr std::pair<LongScalar, LongScalar> estimateIntegralImpl(const Function& f, const Scalar xmin, const Scalar xmax);
	
	inline constexpr Scalar getMaxDeltaXImpl(const Scalar xmin, const Scalar xmax) const { return (xmax - xmin)*misc::maxDiff(std::span{s_xi}); } 
private:
	std::array<LongScalar, 15> m_fx15;
	
	static constexpr inline std::array<Scalar, 15> s_xi{
		Scalar(-0.98799251802048542849), -Scalar(0.93727339240070590431), -Scalar(0.84820658341042721620), Scalar(-0.72441773136017004742), Scalar(-0.57097217260853884754), 
		Scalar(-0.39415134707756336990), -Scalar(0.20119409399743452230),  Scalar(                      ), Scalar( 0.20119409399743452230), Scalar( 0.39415134707756336990), 
		Scalar( 0.57097217260853884754),  Scalar(0.72441773136017004742),  Scalar(0.84820658341042721620), Scalar( 0.93727339240070590431), Scalar( 0.98799251802048542849)};
	
	static constexpr inline std::array<Scalar, 15> s_wi15{
		Scalar(0.030753241996117268355), Scalar(0.070366047488108124709), Scalar(0.10715922046717193501), Scalar(0.13957067792615431445 ), Scalar(0.16626920581699393355 ), 
		Scalar(0.18616100001556221103 ), Scalar(0.19843148532711157646 ), Scalar(0.20257824192556127288), Scalar(0.19843148532711157646 ), Scalar(0.18616100001556221103 ), 
		Scalar(0.16626920581699393355 ), Scalar(0.13957067792615431445 ), Scalar(0.10715922046717193501), Scalar(0.070366047488108124709), Scalar(0.030753241996117268355)};
		
	static constexpr inline std::array<Scalar, 15> s_wi14{
		Scalar(0.042948056434679514011 ), Scalar(0.028746310200837528204), Scalar(0.18519843621047418522), Scalar(0.023655483418631419965), Scalar(0.31694007279358916956  ), 
		Scalar(0.0076858378839750032222), Scalar(0.39482580305781317982 ), Scalar(                      ), Scalar(0.39482580305781317982 ), Scalar(0.0076858378839750032222), 
		Scalar(0.31694007279358916956  ), Scalar(0.023655483418631419965), Scalar(0.18519843621047418522), Scalar(0.028746310200837528204), Scalar(0.042948056434679514011 )};
		
	static constexpr inline std::array<Scalar, 15> s_wi06{
		Scalar(                      ), Scalar(0.21431521897243154264 ), Scalar{}, Scalar(0.062261803859627636066), Scalar{}                     , 
		Scalar(0.72342297716794082129), Scalar(                       ), Scalar{}, Scalar(                       ), Scalar(0.72342297716794082129), 
		Scalar(                      ), Scalar(0.062261803859627636066), Scalar{}, Scalar(0.21431521897243154264 ), Scalar{}};	
};

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature_impl.hpp>

#endif // LNIT_GAUSS_LEGENDRE_ADAPTIVE_QUADRATURE_HPP
