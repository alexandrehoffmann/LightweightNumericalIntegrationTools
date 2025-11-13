#ifndef LNIT_MISC_FP_COMPARATOR_HPP
#define LNIT_MISC_FP_COMPARATOR_HPP

#include <cmath>

namespace LNIT
{
namespace misc
{

template<class Scalar>
class FPComparator
{
public:
	FPComparator(const Scalar& epsilon) : m_epsilon(epsilon) {}
	
	constexpr bool approxEq              (const Scalar& a, const Scalar& b) const { return std::abs(a - b) <= std::max(std::abs(a), std::abs(b))*m_epsilon; }
	constexpr bool essentEq              (const Scalar& a, const Scalar& b) const { return std::abs(a - b) <= std::min(std::abs(a), std::abs(b))*m_epsilon; }
	constexpr bool definitelyGreaterThan (const Scalar& a, const Scalar& b) const { return         (a - b) >  std::max(std::abs(a), std::abs(b))*m_epsilon; }
	constexpr bool definitelyLessThan    (const Scalar& a, const Scalar& b) const { return         (b - a) >  std::max(std::abs(a), std::abs(b))*m_epsilon; }
private:
	Scalar m_epsilon;
};

} // namespace misc
} // namespace LNIT

#endif // LNIT_MISC_FP_COMPARATOR_HPP
