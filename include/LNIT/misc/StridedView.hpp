#ifndef LNIT_MISC_STRIDED_VIEW_HPP
#define LNIT_MISC_STRIDED_VIEW_HPP

#include <LNIT/misc/StridedIterator.hpp>

namespace LNIT
{
namespace misc
{
	
template<typename Values, typename Stride>
class StridedView
{
	using ValuesIterator       = typename Values::iterator;
	using const_ValuesIterator = typename Values::const_iterator;
		
    using Iterator       = StridedIterator<ValuesIterator, Stride>;
    using const_Iterator = StridedIterator<const_ValuesIterator, Stride>;
    
    using ReverseIterator       = StridedIterator<ValuesIterator, Stride>;
    using const_ReverseIterator = StridedIterator<const_ValuesIterator, Stride>;
public:
    StridedView(Values& values, const Stride stride) : m_values(values), m_stride(stride) {}
    
	Iterator begin() { return Iterator(m_values.begin(), m_stride); }
	Iterator end()   { return Iterator(m_values.begin(), m_stride); }

	const_Iterator begin() const { return const_Iterator(m_values.begin(), m_stride); }
	const_Iterator end()   const { return const_Iterator(m_values.begin(), m_stride); }

	const_Iterator cbegin() const { return const_Iterator(m_values.begin(), m_stride); }
	const_Iterator cend()   const { return const_Iterator(m_values.begin(), m_stride); }

	ReverseIterator rbegin() { return ReverseIterator(m_values.begin(), m_stride); }
	ReverseIterator rend()   { return ReverseIterator(m_values.begin(), m_stride); }

	const_ReverseIterator rbegin() const { return const_ReverseIterator(m_values.begin(), m_stride); }
	const_ReverseIterator rend()   const { return const_ReverseIterator(m_values.begin(), m_stride); }

	const_ReverseIterator crbegin() const { return const_ReverseIterator(m_values.begin(), m_stride); }
	const_ReverseIterator crend()   const { return const_ReverseIterator(m_values.begin(), m_stride); }
private:
	Values& m_values;
	Stride  m_stride;
};
	
} // namespace misc
} // namespace LNIT

#endif // LNIT_MISC_STRIDED_VIEW_HPP
