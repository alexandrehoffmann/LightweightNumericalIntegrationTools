#ifndef LNIT_MISC_NUMERIC_IMPL_HPP
#define LNIT_MISC_NUMERIC_IMPL_HPP

#include <span>
#include <cstddef>

#include <LNIT/misc/Numeric.hpp>

namespace LNIT
{
namespace misc
{

template<typename T, size_t N> 
constexpr T maxDiff(std::span<T, N> x)
{
	using Mutable = std::remove_cv_t<T>;
	
	Mutable res = std::abs(x[1] - x[0]);
	for (size_t i=1; i!=size_t(x.size()-1); ++i)
	{
		res = std::max(res, std::abs(x[i+1] - x[i]));
	}
	
	return res;
}

} // namespace misc
} // namespace LNIT

#endif // LNIT_MISC_NUMERIC_IMPL_HPP
