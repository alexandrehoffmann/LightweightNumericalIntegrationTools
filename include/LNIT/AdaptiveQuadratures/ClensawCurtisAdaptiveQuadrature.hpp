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

	static constexpr std::array<Scalar, 9>  s_wi09 = {
		1.0/72.0,
		65.0/252.0 - 8.0/105.0*M_SQRT2,
		347.0/1260.0,
		(8.0/105.0)*M_SQRT2 + 65.0/252.0,
		491.0/1260.0,
		(8.0/105.0)*M_SQRT2 + 65.0/252.0,
		347.0/1260.0,
		65.0/252.0 - 8.0/105.0*M_SQRT2,
		1.0/72.0};
	static constexpr std::array<Scalar, 17> s_wi17 = {
		1.0/272.0,
		-16.0/195.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 16.0/2145.0*M_SQRT2 - 16.0/3465.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 257.0/2040.0,
		5483.0/42840.0 - 1688.0/45045.0*M_SQRT2,
		-16.0/195.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (16.0/3465.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (16.0/2145.0)*M_SQRT2 + 257.0/2040.0,
		848261.0/6126120.0,
		-16.0/3465.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (16.0/2145.0)*M_SQRT2 + (16.0/195.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 257.0/2040.0,
		(1688.0/45045.0)*M_SQRT2 + 5483.0/42840.0,
		-16.0/2145.0*M_SQRT2 + (16.0/3465.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (16.0/195.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + 257.0/2040.0,
		1200229.0/6126120.0,
		-16.0/2145.0*M_SQRT2 + (16.0/3465.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (16.0/195.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + 257.0/2040.0,
		(1688.0/45045.0)*M_SQRT2 + 5483.0/42840.0,
		-16.0/3465.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (16.0/2145.0)*M_SQRT2 + (16.0/195.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 257.0/2040.0,
		848261.0/6126120.0,
		-16.0/195.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (16.0/3465.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (16.0/2145.0)*M_SQRT2 + 257.0/2040.0,
		5483.0/42840.0 - 1688.0/45045.0*M_SQRT2,
		-16.0/195.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 16.0/2145.0*M_SQRT2 - 16.0/3465.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 257.0/2040.0,
		1.0/272.0};
	static constexpr std::array<Scalar, 33> s_wi33 = {
		1.0/1056.0,
		-112.0/2697.0*std::cos((1.0/16.0)*M_PI) - 32.0/3915.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 16.0/4725.0*std::cos((3.0/16.0)*M_PI) - 32.0/36225.0*M_SQRT2 - 16.0/15939.0*std::cos((5.0/16.0)*M_PI) - 32.0/57057.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 16.0/62985.0*std::cos((7.0/16.0)*M_PI) + 1025.0/16368.0,
		-2308976.0/56623515.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 271984.0/74459385.0*M_SQRT2 - 2672.0/1195425.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 87467.0/1391280.0,
		-112.0/2697.0*std::cos((3.0/16.0)*M_PI) - 32.0/3915.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 1.0/3864.0*std::cos((33.0/16.0)*M_PI) - 1.0/5400.0*std::cos((39.0/16.0)*M_PI) + (16.0/62985.0)*std::cos((5.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (1.0/280.0)*std::cos((7.0/16.0)*M_PI) + (1.0/792.0)*std::cos((1.0/16.0)*M_PI) + (32.0/36225.0)*M_SQRT2 + 1025.0/16368.0,
		215335717.0/3359941200.0 - 84737118104.0/4512611027925.0*M_SQRT2,
		-112.0/2697.0*std::cos((5.0/16.0)*M_PI) - 1.0/1560.0*std::cos((35.0/16.0)*M_PI) - 32.0/57057.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 16.0/15939.0*std::cos((7.0/16.0)*M_PI) - 1.0/5400.0*std::cos((65.0/16.0)*M_PI) + (1.0/2584.0)*std::cos((3.0/16.0)*M_PI) + (32.0/36225.0)*M_SQRT2 + (32.0/3915.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/280.0)*std::cos((1.0/16.0)*M_PI) + 1025.0/16368.0,
		-2308976.0/56623515.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (2672.0/1195425.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (271984.0/74459385.0)*M_SQRT2 + 87467.0/1391280.0,
		-112.0/2697.0*std::cos((7.0/16.0)*M_PI) - 32.0/36225.0*M_SQRT2 - 1.0/792.0*std::cos((35.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/3864.0)*std::cos((3.0/16.0)*M_PI) + (16.0/62985.0)*std::cos((1.0/16.0)*M_PI) + (16.0/4725.0)*std::cos((5.0/16.0)*M_PI) + (32.0/3915.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + 1025.0/16368.0,
		5008383140437.0/72201776446800.0,
		-16.0/4725.0*std::cos((5.0/16.0)*M_PI) - 32.0/36225.0*M_SQRT2 - 16.0/62985.0*std::cos((1.0/16.0)*M_PI) - 1.0/3864.0*std::cos((99.0/16.0)*M_PI) - 1.0/7192.0*std::cos((135.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/792.0)*std::cos((3.0/16.0)*M_PI) + (32.0/3915.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (1.0/24.0)*std::cos((7.0/16.0)*M_PI) + 1025.0/16368.0,
		-2672.0/1195425.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (271984.0/74459385.0)*M_SQRT2 + (2308976.0/56623515.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 87467.0/1391280.0,
		-1.0/280.0*std::cos((33.0/16.0)*M_PI) - 32.0/57057.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 1.0/2584.0*std::cos((99.0/16.0)*M_PI) - 1.0/7192.0*std::cos((165.0/16.0)*M_PI) + (1.0/5400.0)*std::cos((1.0/16.0)*M_PI) + (16.0/15939.0)*std::cos((7.0/16.0)*M_PI) + (1.0/1560.0)*std::cos((3.0/16.0)*M_PI) + (32.0/36225.0)*M_SQRT2 + (32.0/3915.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/24.0)*std::cos((5.0/16.0)*M_PI) + 1025.0/16368.0,
		(84737118104.0/4512611027925.0)*M_SQRT2 + 215335717.0/3359941200.0,
		-32.0/3915.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 1.0/792.0*std::cos((65.0/16.0)*M_PI) - 1.0/280.0*std::cos((39.0/16.0)*M_PI) - 16.0/62985.0*std::cos((5.0/16.0)*M_PI) - 1.0/7192.0*std::cos((195.0/16.0)*M_PI) + (1.0/5400.0)*std::cos((7.0/16.0)*M_PI) + (1.0/3864.0)*std::cos((1.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (32.0/36225.0)*M_SQRT2 + (1.0/24.0)*std::cos((3.0/16.0)*M_PI) + 1025.0/16368.0,
		-271984.0/74459385.0*M_SQRT2 + (2672.0/1195425.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (2308976.0/56623515.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + 87467.0/1391280.0,
		-32.0/3915.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 32.0/36225.0*M_SQRT2 - 32.0/57057.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 1.0/5400.0*std::cos((195.0/16.0)*M_PI) - 1.0/3864.0*std::cos((165.0/16.0)*M_PI) - 1.0/7192.0*std::cos((225.0/16.0)*M_PI) - 1.0/2584.0*std::cos((135.0/16.0)*M_PI) + (1.0/1560.0)*std::cos((7.0/16.0)*M_PI) + (1.0/792.0)*std::cos((5.0/16.0)*M_PI) + (1.0/280.0)*std::cos((3.0/16.0)*M_PI) + (1.0/24.0)*std::cos((1.0/16.0)*M_PI) + 1025.0/16368.0,
		7084256570389.0/72201776446800.0,
		-32.0/3915.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 32.0/36225.0*M_SQRT2 - 32.0/57057.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 1.0/5400.0*std::cos((195.0/16.0)*M_PI) - 1.0/3864.0*std::cos((165.0/16.0)*M_PI) - 1.0/7192.0*std::cos((225.0/16.0)*M_PI) - 1.0/2584.0*std::cos((135.0/16.0)*M_PI) + (1.0/1560.0)*std::cos((7.0/16.0)*M_PI) + (1.0/792.0)*std::cos((5.0/16.0)*M_PI) + (1.0/280.0)*std::cos((3.0/16.0)*M_PI) + (1.0/24.0)*std::cos((1.0/16.0)*M_PI) + 1025.0/16368.0,
		-271984.0/74459385.0*M_SQRT2 + (2672.0/1195425.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (2308976.0/56623515.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + 87467.0/1391280.0,
		-32.0/3915.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 1.0/792.0*std::cos((65.0/16.0)*M_PI) - 1.0/280.0*std::cos((39.0/16.0)*M_PI) - 16.0/62985.0*std::cos((5.0/16.0)*M_PI) - 1.0/7192.0*std::cos((195.0/16.0)*M_PI) + (1.0/5400.0)*std::cos((7.0/16.0)*M_PI) + (1.0/3864.0)*std::cos((1.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (32.0/36225.0)*M_SQRT2 + (1.0/24.0)*std::cos((3.0/16.0)*M_PI) + 1025.0/16368.0,
		(84737118104.0/4512611027925.0)*M_SQRT2 + 215335717.0/3359941200.0,
		-1.0/280.0*std::cos((33.0/16.0)*M_PI) - 32.0/57057.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 1.0/2584.0*std::cos((99.0/16.0)*M_PI) - 1.0/7192.0*std::cos((165.0/16.0)*M_PI) + (1.0/5400.0)*std::cos((1.0/16.0)*M_PI) + (16.0/15939.0)*std::cos((7.0/16.0)*M_PI) + (1.0/1560.0)*std::cos((3.0/16.0)*M_PI) + (32.0/36225.0)*M_SQRT2 + (32.0/3915.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/24.0)*std::cos((5.0/16.0)*M_PI) + 1025.0/16368.0,
		-2672.0/1195425.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (271984.0/74459385.0)*M_SQRT2 + (2308976.0/56623515.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 87467.0/1391280.0,
		-16.0/4725.0*std::cos((5.0/16.0)*M_PI) - 32.0/36225.0*M_SQRT2 - 16.0/62985.0*std::cos((1.0/16.0)*M_PI) - 1.0/3864.0*std::cos((99.0/16.0)*M_PI) - 1.0/7192.0*std::cos((135.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/792.0)*std::cos((3.0/16.0)*M_PI) + (32.0/3915.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (1.0/24.0)*std::cos((7.0/16.0)*M_PI) + 1025.0/16368.0,
		5008383140437.0/72201776446800.0,
		-112.0/2697.0*std::cos((7.0/16.0)*M_PI) - 32.0/36225.0*M_SQRT2 - 1.0/792.0*std::cos((35.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/3864.0)*std::cos((3.0/16.0)*M_PI) + (16.0/62985.0)*std::cos((1.0/16.0)*M_PI) + (16.0/4725.0)*std::cos((5.0/16.0)*M_PI) + (32.0/3915.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + 1025.0/16368.0,
		-2308976.0/56623515.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (2672.0/1195425.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (271984.0/74459385.0)*M_SQRT2 + 87467.0/1391280.0,
		-112.0/2697.0*std::cos((5.0/16.0)*M_PI) - 1.0/1560.0*std::cos((35.0/16.0)*M_PI) - 32.0/57057.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 16.0/15939.0*std::cos((7.0/16.0)*M_PI) - 1.0/5400.0*std::cos((65.0/16.0)*M_PI) + (1.0/2584.0)*std::cos((3.0/16.0)*M_PI) + (32.0/36225.0)*M_SQRT2 + (32.0/3915.0)*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + (1.0/280.0)*std::cos((1.0/16.0)*M_PI) + 1025.0/16368.0,
		215335717.0/3359941200.0 - 84737118104.0/4512611027925.0*M_SQRT2,
		-112.0/2697.0*std::cos((3.0/16.0)*M_PI) - 32.0/3915.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 1.0/3864.0*std::cos((33.0/16.0)*M_PI) - 1.0/5400.0*std::cos((39.0/16.0)*M_PI) + (16.0/62985.0)*std::cos((5.0/16.0)*M_PI) + (32.0/57057.0)*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) + (1.0/280.0)*std::cos((7.0/16.0)*M_PI) + (1.0/792.0)*std::cos((1.0/16.0)*M_PI) + (32.0/36225.0)*M_SQRT2 + 1025.0/16368.0,
		-2308976.0/56623515.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 271984.0/74459385.0*M_SQRT2 - 2672.0/1195425.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) + 87467.0/1391280.0,
		-112.0/2697.0*std::cos((1.0/16.0)*M_PI) - 32.0/3915.0*std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0) - 16.0/4725.0*std::cos((3.0/16.0)*M_PI) - 32.0/36225.0*M_SQRT2 - 16.0/15939.0*std::cos((5.0/16.0)*M_PI) - 32.0/57057.0*std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2) - 16.0/62985.0*std::cos((7.0/16.0)*M_PI) + 1025.0/16368.0,
		1.0/1056.0};
	static constexpr std::array<Scalar, 33> s_xi = {
		 1,
		 std::cos((1.0/32.0)*M_PI),
		 std::cos((1.0/16.0)*M_PI),
		 std::cos((3.0/32.0)*M_PI),
		 std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0),
		 std::cos((5.0/32.0)*M_PI),
		 std::cos((3.0/16.0)*M_PI),
		 std::cos((7.0/32.0)*M_PI),
		 M_SQRT1_2,
		 std::cos((9.0/32.0)*M_PI),
		 std::cos((5.0/16.0)*M_PI),
		 std::cos((11.0/32.0)*M_PI),
		 std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2),
		 std::cos((13.0/32.0)*M_PI),
		 std::cos((7.0/16.0)*M_PI),
		 std::cos((15.0/32.0)*M_PI),
		 0,
		-std::cos((15.0/32.0)*M_PI),
		-std::cos((7.0/16.0)*M_PI),
		-std::cos((13.0/32.0)*M_PI),
		-std::sqrt(1.0/2.0 - 1.0/4.0*M_SQRT2),
		-std::cos((11.0/32.0)*M_PI),
		-std::cos((5.0/16.0)*M_PI),
		-std::cos((9.0/32.0)*M_PI),
		-M_SQRT1_2,
		-std::cos((7.0/32.0)*M_PI),
		-std::cos((3.0/16.0)*M_PI),
		-std::cos((5.0/32.0)*M_PI),
		-std::sqrt((1.0/4.0)*M_SQRT2 + 1.0/2.0),
		-std::cos((3.0/32.0)*M_PI),
		-std::cos((1.0/16.0)*M_PI),
		-std::cos((1.0/32.0)*M_PI),
		-1};
};
	
} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/ClensawCurtisAdaptiveQuadrature_impl.hpp>

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_HPP
