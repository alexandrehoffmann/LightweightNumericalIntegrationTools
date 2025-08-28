#ifndef LNIT_MISC_ARRAY_VIEW_HPP
#define LNIT_MISC_ARRAY_VIEW_HPP

#include <LNIT/misc/ArrayViewIterator.hpp>

namespace LNIT
{
namespace misc
{

/**
 * @brief Provides a view into a value container through an index container.
 *
 * This class acts as an adapter that allows iteration over a subset or reordered
 * sequence of values, as defined by an index container. The indices determine
 * which elements from the underlying value container are exposed and in what order.
 *
 * @tparam Values Container type holding the underlying values (e.g. std::vector<T>).
 * @tparam Idxs   Container or Range type holding indices (e.g. std::vector<std::size_t>, std::ranges::iota_view, etc.).
 */
template<typename Values, typename Idxs>
class ArrayView
{
	using ValuesIterator       = typename Values::iterator;
	using const_ValuesIterator = typename Values::const_iterator;
	using IdxsIterator         = typename Idxs::const_iterator;
    using ReverseIdxsIterator  = typename Idxs::const_reverse_iterator;

    using Iterator       = ArrayViewIterator<ValuesIterator, IdxsIterator>;
    using const_Iterator = ArrayViewIterator<const_ValuesIterator, IdxsIterator>;
    
    using ReverseIterator       = ArrayViewIterator<ValuesIterator, ReverseIdxsIterator>;
    using const_ReverseIterator = ArrayViewIterator<const_ValuesIterator, ReverseIdxsIterator>;
public:
	/**
     * @brief Constructs an ArrayView.
     *
     * @param values Reference to the value container.
     * @param idxs   Reference to the index container.
     */
	ArrayView(Values& values, Idxs& idxs) : m_values(values), m_idxs(idxs) {}
	
	Iterator begin() { return Iterator(m_values.begin(), m_idxs.begin()); }
    Iterator end()   { return Iterator(m_values.begin(), m_idxs.end()); }

    const_Iterator begin() const { return const_Iterator(m_values.begin(), m_idxs.begin()); }
    const_Iterator end()   const { return const_Iterator(m_values.begin(), m_idxs.end()); }

    const_Iterator cbegin() const { return const_Iterator(m_values.begin(), m_idxs.begin()); }
    const_Iterator cend()   const { return const_Iterator(m_values.begin(), m_idxs.end()); }

	ReverseIterator rbegin() { return ReverseIterator(m_values.begin(), m_idxs.rbegin()); }
    ReverseIterator rend()   { return ReverseIterator(m_values.begin(), m_idxs.rend()); }

    const_ReverseIterator rbegin() const { return const_ReverseIterator(m_values.begin(), m_idxs.crbegin()); }
    const_ReverseIterator rend()   const { return const_ReverseIterator(m_values.begin(), m_idxs.crend()); }

    const_ReverseIterator crbegin() const { return const_ReverseIterator(m_values.begin(), m_idxs.crbegin()); }
    const_ReverseIterator crend()   const { return const_ReverseIterator(m_values.begin(), m_idxs.crend()); }
private:
	Values& m_values;
	Idxs&   m_idxs;
};
	
} // namespace misc
} // namespace LNIT

#endif // LNIT_MISC_ARRAY_VIEW_HPP
