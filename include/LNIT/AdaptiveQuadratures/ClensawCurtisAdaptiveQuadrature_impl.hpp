#ifndef LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
#define LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP

#include <LNIT/AdaptiveQuadratures/ClensawCurtisAdaptiveQuadrature.hpp>

#include <LNIT/misc/StridedIterator.hpp>

namespace LNIT
{

//// ndoes and weights initialization ////

template<typename T, typename TT> 
const std::array<typename ClensawCurtisAdaptiveQuadrature<T,TT>::Scalar, 9> ClensawCurtisAdaptiveQuadrature<T,TT>::s_wi09 = {
		1.0/72.0,
		65.0/252.0 - 8.0/105.0*M_SQRT2,
		347.0/1260.0,
		(8.0/105.0)*M_SQRT2 + 65.0/252.0,
		491.0/1260.0,
		(8.0/105.0)*M_SQRT2 + 65.0/252.0,
		347.0/1260.0,
		65.0/252.0 - 8.0/105.0*M_SQRT2,
		1.0/72.0};

template<typename T, typename TT> 
const std::array<typename ClensawCurtisAdaptiveQuadrature<T,TT>::Scalar, 17> ClensawCurtisAdaptiveQuadrature<T,TT>::s_wi17 = {
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

template<typename T, typename TT> 	
const std::array<typename ClensawCurtisAdaptiveQuadrature<T,TT>::Scalar, 33> ClensawCurtisAdaptiveQuadrature<T,TT>::s_wi33 = {
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

template<typename T, typename TT> 	
const std::array<typename ClensawCurtisAdaptiveQuadrature<T,TT>::Scalar, 33> ClensawCurtisAdaptiveQuadrature<T,TT>::s_xi = {
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

//// explicit template instanciations ////

extern template class ClensawCurtisAdaptiveQuadrature<double, double>;
extern template class ClensawCurtisAdaptiveQuadrature<double, long double>;

extern template class ClensawCurtisAdaptiveQuadrature<long double, long double>;

//// method implementations ////

template<typename T, typename TT> template<class Function>
auto ClensawCurtisAdaptiveQuadrature<T,TT>::estimateIntegral_impl(const Function& f, const Scalar xmin, const Scalar xmax) -> std::pair<LongScalar, LongScalar>
{	
    for (Size i=0; i!=s_xi.size(); ++i)
    {
        const Scalar x = 0.5*(s_xi[i]*(xmax - xmin) + (xmax + xmin));
        m_fx[i] = f(x);
    }

    const LongScalar I33 = 0.5*(xmax - xmin)*std::inner_product(s_wi33.cbegin(), s_wi33.cend(), m_fx.cbegin(), LongScalar(0));
    const LongScalar I17 = 0.5*(xmax - xmin)*std::inner_product(s_wi17.cbegin(), s_wi17.cend(), misc::StridedIterator(m_fx.cbegin(), 2), LongScalar(0));
    const LongScalar I09 = 0.5*(xmax - xmin)*std::inner_product(s_wi09.cbegin(), s_wi09.cend(), misc::StridedIterator(m_fx.cbegin(), 4), LongScalar(0));

    const LongScalar err1 = std::abs(I33 - I17);
	const LongScalar err2 = std::abs(I33 - I09);

    return std::make_pair(I33, err2 == 0 ? LongScalar(0) : err2*(err1 / err2)*(err1 / err2));
}

} // namespace LNIT

#endif // LNIT_CLENSAW_CURTIS_ADAPTIVE_QUADRATURE_IMPL_HPP
