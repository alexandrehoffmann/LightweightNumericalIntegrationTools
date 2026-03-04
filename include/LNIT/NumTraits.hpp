#ifndef LNIT_NUM_TRAITS_HPP
#define LNIT_NUM_TRAITS_HPP

#include <concepts>
#include <limits>

namespace LNIT
{

template<typename T> struct NumTraits;

template<std::floating_point T>
struct NumTraits<T>
{
	static constexpr T epsilon = std::numeric_limits<T>::epsilon();
	static constexpr T NaN     = std::numeric_limits<T>::quiet_NaN();
};

} // namespace LNIT

#endif // LNIT_NUM_TRAITS_HPP
