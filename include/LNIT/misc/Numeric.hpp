#ifndef LNIT_MISC_NUMERIC_HPP
#define LNIT_MISC_NUMERIC_HPP

#include <span>
#include <cstddef>

namespace LNIT
{
namespace misc
{

template<typename T, size_t N> constexpr T maxDiff(std::span<T, N> x);

} // namespace misc
} // namespace LNIT

#include <LNIT/misc/Numeric_impl.hpp>

#endif // LNIT_MISC_NUMERIC_HPP
