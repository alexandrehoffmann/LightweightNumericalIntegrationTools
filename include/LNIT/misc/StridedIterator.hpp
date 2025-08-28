#ifndef LNIT_MISC_STRIDED_ITERATOR_HPP
#define LNIT_MISC_STRIDED_ITERATOR_HPP

#include <iterator>

namespace LNIT
{
namespace misc
{

/**
 * @brief Iterator adaptor that advances an underlying iterator by a fixed stride.
 *
 * Useful for evaluating quadrature rules with reduced nodes (e.g., skipping points).
 *
 * @todo remove when moving to c++23
 *
 * @tparam Iterator The underlying iterator type.
 * @tparam Stride The stride type (e.g., int).
 */
template<class Iterator, typename Stride>
class StridedIterator
{
    using iterator_category = std::input_iterator_tag;
    using difference_type   = typename std::iterator_traits<Iterator>::difference_type;
    using value_type        = typename std::iterator_traits<Iterator>::value_type;
    using pointer           = typename std::iterator_traits<Iterator>::pointer;
    using reference         = typename std::iterator_traits<Iterator>::reference;
public:
    StridedIterator(Iterator firstValue, Stride stride) : m_firstValue(firstValue), m_stride(stride) {}

    reference operator*()  const { return *m_firstValue; } 
	pointer   operator->()       { return std::addressof(*m_firstValue); }

    /// @brief Prefix increment.
	StridedIterator& operator++() { m_firstValue += m_stride; return *this; }

	/// @brief Prefix increment.
	StridedIterator operator++(int) { StridedIterator tmp = *this; ++(*this); return tmp; }

	friend bool operator== (const StridedIterator& a, const StridedIterator& b) { return a.m_firstValue == b.m_firstValue; }
	friend bool operator!= (const StridedIterator& a, const StridedIterator& b) { return a.m_firstValue != b.m_firstValue; }  
private:
    Iterator m_firstValue;
    Stride   m_stride;
};
	
} // namespace misc
} // namespace LNIT

#endif // LNIT_MISC_STRIDED_ITERATOR_HPP
