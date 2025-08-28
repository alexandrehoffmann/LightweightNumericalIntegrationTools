#ifndef LNIT_ADAPTIVE_QUADRATURE_BASE_HPP
#define LNIT_ADAPTIVE_QUADRATURE_BASE_HPP

#include <cstdio>
#include <array>
#include <vector>

namespace LNIT
{

template<class Derived> class AdaptiveQuadratureBase;

template<class Derived> struct AdaptiveQuadratureTraits;

template<class Derived>
class AdaptiveQuadratureBase
{
	using DerivedTraits = AdaptiveQuadratureTraits<Derived>;
public:
	using Size       = typename DerivedTraits::Size;       ///<  @brief Type for iteration counters.
	using Scalar     = typename DerivedTraits::Scalar;     ///<  @brief Floating point type for integration (e.g., double).
	using LongScalar = typename DerivedTraits::LongScalar; ///<  @brief Higher precision type for accumulation (e.g., long double).
	using Interval = std::pair<Scalar, Scalar>;            ///<  @brief Type representing an integration interval.
	
	      Derived& derived_cast()       { return static_cast<      Derived&>(*this); }
	const Derived& derived_cast() const { return static_cast<const Derived&>(*this); }
	
	/**
	 * @brief Construct quadrature object.
	 * @param maxIt Maximum number of iterations.
	 * @param tol Relative tolerance for convergence.
	 */
	AdaptiveQuadratureBase(const Size maxIt = 20000, const Scalar tol = std::numeric_limits<Scalar>::epsilon());

	/**
	 * @brief Estimate integral and error on [xmin, xmax].
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound.
	 * @param xmax Upper bound.
	 * 
	 * Delegates to Derived::estimateIntegral_impl.
	 * 
	 * @return Pair (integral, estimated error).
	 */
	template<class Function> std::pair<LongScalar, LongScalar> estimateIntegral(const Function& f, const Scalar xmin, const Scalar xmax) { return derived_cast().estimateIntegral_impl(f, xmin, xmax); }

	/**
	 * @brief Perform adaptive quadrature on [xmin, xmax].
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound of interval.
	 * @param xmax Upper bound of interval.
	 * @return Approximation of the integral.
	 */
	template<class Function> Scalar integrate(const Function& f, const Scalar xmin, const Scalar xmax);
	
	/**
	 * @brief Perform adaptive quadrature on (-inf, xmax].
	 * 
	 * First tries to find xmin such that:
	 * \f[
	 * \int_{-\infty}^{xmin} f(x) dx \approx 0.
	 * \f]
	 * And then addapt the quadrature over [xmin, xmax].
	 * 
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound of interval.
	 * @param xmax Upper bound of interval.
	 * @return Approximation of the integral.
	 */
	template<class Function> Scalar integrateLeftInfinite(const Function& f, const Scalar xmax);
	
	/**
	 * @brief Perform adaptive quadrature on [xmin, inf).
	 * 
	 * First tries to find xmax such that:
	 * \f[
	 * \int_{xmax}^{\infty} f(x) dx \approx 0.
	 * \f]
	 * And then addapt the quadrature over [xmin, xmax].
	 * 
	 * @tparam Function Callable with signature Scalar f(Scalar).
	 * @param f Function to integrate.
	 * @param xmin Lower bound of interval.
	 * @param xmax Upper bound of interval.
	 * @return Approximation of the integral.
	 */
	template<class Function> Scalar integrateRightInfinite(const Function& f, const Scalar xmin);
	
	/**
     * @brief Perform adaptive quadrature on (-inf, inf).
     * 
     * First tries to find xmin such that:
	 * \f[
	 * \int_{-\infty}^{xmin} f(x) dx \approx 0.
	 * \f]
	 * Then tries to find xmax such that:
	 * \f[
	 * \int_{xmax}^{\infty} f(x) dx \approx 0.
	 * \f]
	 * Finally addapt the quadrature over [xmin, xmax].
     */
	template<class Function> Scalar integrate(const Function& f);
	
    Size   getMaxIt() const { return m_maxIt; } ///<  @brief Maximum iterations allowed.
	Size   getNits()  const { return m_it; }    ///<  @brief Number of iterations performed.
	Scalar getTol()   const { return m_tol; }   ///<  @brief Tolerance of the quadrature.
	
	bool hasConverged() const { return m_hasConverged; } ///<  @brief Whether convergence was achieved.
	
	void setMaxIt(const Size maxIt) { m_maxIt = maxIt; } ///<  @brief Set maximum number of iterations.
	void setTol(const Scalar& tol)  { m_tol   = tol; }   ///<  @brief Set the tolerance of the quadrature.
	
	void setOutput(std::FILE* out) { m_out = out; } ///<  @brief Redirect output to a file for logging progress.
private:
	std::vector<Interval>   m_intervals;
	std::vector<LongScalar> m_subIntergrals;
	std::vector<LongScalar> m_subIntergralsErr;

	Size   m_maxIt;
	Size   m_it;
	Scalar m_tol;
	bool   m_hasConverged;
	
	std::FILE* m_out = nullptr;
};

template<class T> struct IsAdaptiveQuadrature : std::bool_constant< std::is_base_of<AdaptiveQuadratureBase<T>, T>::value > {};  ///<  @brief Trait to determine if a type derives from AdaptiveQuadratureBase.

template<class T> concept IsAdaptiveQuadrature_concept = IsAdaptiveQuadrature<T>::value;

} // namespace LNIT

#include <LNIT/AdaptiveQuadratures/AdaptiveQuadratureBase_impl.hpp>

#endif // LNIT_ADAPTIVE_QUADRATURE_BASE_HPP
