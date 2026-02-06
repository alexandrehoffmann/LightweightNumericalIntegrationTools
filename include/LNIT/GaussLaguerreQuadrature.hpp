#ifndef LNIT_GAUSS_LAGUERRE_QUADRATURE_HPP
#define LNIT_GAUSS_LAGUERRE_QUADRATURE_HPP

#include <array>

namespace LNIT
{

/**
 * @class GaussLaguerreQuadrature
 * @brief Template class implementing Gauss-Laguerre quadrature.
 *
 * @tparam Scalar Numeric type for computation (e.g., float, double).
 * @tparam LongScalar Optional extended precision type for accumulation (defaults to Scalar).
 *
 * This class provides methods to numerically approximate integrals of functions defined
 * on semi-infinite intervals using Gauss-Laguerre quadrature rules.
 */
template<typename Scalar, typename LongScalar=Scalar>
class GaussLaguerreQuadrature
{
public:
	using Size = unsigned int; ///< Type alias for sizes/indices.

	GaussLaguerreQuadrature() {} ///<  @brief Default constructor.

	/**
	 * @brief Approximate integral over the left semi-infinite interval.
	 *
	 * Computes an approximation of the integral:
	 * \f[
	 * \int_{-\infty}^{a} f(x) dx
	 * \f]
	 * using Gauss-Laguerre quadrature.
	 *
	 * @tparam Function Type of the callable object (e.g., lambda, functor).
	 * @param f Integrand function \f$ f(x) \f$.
	 * @param a Upper bound of integration (default = 0).
	 * @return Approximated integral value.
	 */
	template<class Function> Scalar integrateLeftInfinite  (const Function& f, const Scalar a = 0.) const;
	/**
	 * @brief Approximate integral over the right semi-infinite interval.
	 *
	 * Computes an approximation of the integral:
	 * \f[
	 * \int_{a}^{+\infty} f(x) dx
	 * \f]
	 * using Gauss-Laguerre quadrature.
	 *
	 * @tparam Function Type of the callable object (e.g., lambda, functor).
	 * @param f Integrand function \f$ f(x) \f$.
	 * @param a Lower bound of integration (default = 0).
	 * @return Approximated integral value.
	 */
	template<class Function> Scalar integrateRightInfinite (const Function& f, const Scalar a = 0.) const;
private:
	static constexpr std::array<Scalar, 33> s_wi = {
		Scalar(0.11077730587320757274), Scalar(0.25810528128189475158), Scalar(0.40622176868437369247), Scalar(0.55526230959922306292), Scalar(0.70555738765958285661), 
		Scalar(0.85746574519748546260), Scalar(1.0113689523807972871 ), Scalar(1.1676756607019108316 ), Scalar(1.3268279171086353030 ), Scalar(1.4893088816572998965 ), 
		Scalar(1.6556521352112563334 ), Scalar(1.8264530047051926656 ), Scalar(2.0023825182904239873 ), Scalar(2.1842048313456294514 ), Scalar(2.3727992840911420803 ), 
		Scalar(2.5691887155313371033 ), Scalar(2.7745763467843133346 ), Scalar(2.9903945894647637769 ), Scalar(3.2183707497169848072 ), Scalar(3.4606171628677355124 ), 
		Scalar(3.7197574810055936975 ), Scalar(3.9991078929934509694 ), Scalar(4.3029443879423609989 ), Scalar(4.6369096377011849903 ), Scalar(5.0086560370564072925 ), 
		Scalar(5.4289084747474177774 ), Scalar(5.9133193778216559086 ), Scalar(6.4859354139778471922 ), Scalar(7.1862724728538016301 ), Scalar(8.0855716526458441488 ), 
		9.3310592514615706576,  11.303622656587531387,  15.566291494163869906};
	static constexpr std::array<Scalar, 33> s_xi = {
		0.043161135617326892192, 0.22751780280337112385, 0.55961665585153988759, 1.0402685077510020538, 1.6705591960757151909, 
		2.4519207958976305465,   3.3861553375880048323,  4.4754594983997714570,  5.7224547202721035227, 7.1302243444001080163, 
		8.7023592306214062489,   10.443013650205982427,  12.356973759350285962,  14.449741681585540238, 16.727639218638322353, 
		19.197936587212446637,   21.869013524928189871,  24.750562906157795643,  27.853851111413356780, 31.192055545575129868, 
		34.780709153538337700,   38.638296717774030225,  42.787072078253479488,  47.254206602993265817, 52.073451901514220267, 
		57.287634541092940075,   62.952565946906630207,  69.143513380109892446,  75.966687014247062344, 83.581637223270880761, 
		92.251139444135101234,   102.47784433682332258,  115.55475644899580731};
};

} // namespace LNIT

#include <LNIT/GaussLaguerreQuadrature_impl.hpp>

#endif // LNIT_GAUSS_LAGUERRE_QUADRATURE_HPP
