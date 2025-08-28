#ifndef LNIT_MISC_ARRAY_VIEW_ITERATOR_HPP
#define LNIT_MISC_ARRAY_VIEW_ITERATOR_HPP

#include <iterator>

namespace LNIT
{
namespace misc
{

/**
 * @brief An iterator that provides indexed access into a value array via an index array.
 *
 * This iterator adapts two iterators: one pointing to the beginning of a value sequence,
 * and one pointing into a sequence of indices. Dereferencing the iterator yields an element
 * from the value sequence at the position specified by the current index.
 *
 * @tparam ValueIterator Iterator type pointing to the underlying value array.
 * @tparam IdxIterator   Iterator type pointing to the index array.
 */
template<typename ValueIterator, typename IdxIterator>
class ArrayViewIterator
{
public:
	using iterator_category = std::bidirectional_iterator_tag;
    using difference_type   = typename std::iterator_traits<ValueIterator>::difference_type;
    using value_type        = typename std::iterator_traits<ValueIterator>::value_type;
    using pointer           = typename std::iterator_traits<ValueIterator>::pointer;
    using reference         = typename std::iterator_traits<ValueIterator>::reference;

	/**
     * @brief Constructs an ArrayViewIterator.
     *
     * @param firstValue Iterator to the beginning of the value array.
     * @param firstIdx   Iterator to the current position in the index array.
     */
	ArrayViewIterator(ValueIterator firstValue, IdxIterator firstIdx) : m_firstValue(firstValue), m_idx(firstIdx) {}
	
	reference operator*() const { return *(m_firstValue + *m_idx); } 
	pointer operator->() { return m_firstValue + *m_idx; }

	ArrayViewIterator& operator++() { m_idx++; return *this; }
    ArrayViewIterator& operator--() { m_idx--; return *this; }  

	ArrayViewIterator operator++(int) { ArrayViewIterator tmp = *this; ++(*this); return tmp; }
    ArrayViewIterator operator--(int) { ArrayViewIterator tmp = *this; --(*this); return tmp; }

	friend bool operator== (const ArrayViewIterator& a, const ArrayViewIterator& b) { return a.m_firstValue == b.m_firstValue and a.m_idx == b.m_idx; }
	friend bool operator!= (const ArrayViewIterator& a, const ArrayViewIterator& b) { return a.m_firstValue != b.m_firstValue or  a.m_idx != b.m_idx; }  
private:
	ValueIterator m_firstValue;
	IdxIterator   m_idx;
};
	
} // namespace misc
} // namespace LNIT

#endif // LNIT_MISC_ARRAY_VIEW_ITERATOR_HPP
