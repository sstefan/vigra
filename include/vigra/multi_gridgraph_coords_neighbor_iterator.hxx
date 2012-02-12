#ifndef VIGRA_MULTI_GRIDGRAPH_COORDS_NEIGHBOR_ITERATOR_HXX
#define VIGRA_MULTI_GRIDGRAPH_COORDS_NEIGHBOR_ITERATOR_HXX

#include "graphs.hxx"

namespace vigra {
    namespace detail {

	template<unsigned int N>
	class CoordsGridGraphNeighborIterator
	{
	public:
	    typedef typename MultiArrayShape<N>::type shape_type;

	    typedef shape_type vertex_descriptor;

	    typedef typename CoupledScanOrderIteratorFactory<N>::coupled_iter_type vertex_iterator;

	    typedef CoordsGridGraphNeighborIterator<N> self_type;

	    typedef vertex_descriptor value_type;

	    typedef value_type*  pointer;
	    typedef value_type const * const_pointer;
	    typedef value_type& reference;
	    typedef value_type const & const_reference;

	    typedef std::ptrdiff_t difference_type;
	    //    typedef std::input_iterator_tag  iterator_category;
	    typedef std::forward_iterator_tag  iterator_category;


	    CoordsGridGraphNeighborIterator() : index_(0),
						neighborhood_(0),
						neighborIndexOffsets_(0),
						neighborIndexLookup_(0)
	    {}

	    CoordsGridGraphNeighborIterator(vertex_descriptor pos,
					    const ArrayVector<shape_type> &neighborhood,
					    const ArrayVector<MultiArrayIndex> &neighborIndexOffsets,
					    const ArrayVector<int> &neighborIndexLookup,
					    bool end=false)
		: pos_(pos),
		  index_(0),
		  neighborhood_(&neighborhood), 
		  neighborIndexOffsets_(&neighborIndexOffsets),
		  neighborIndexLookup_(&neighborIndexLookup)
	    {
		if (end)
		    index_ += (*neighborIndexOffsets_).size();
		else {
		    const shape_type& neighborOffset((*neighborhood_)[index_]);
		    tmp_ = pos_ + neighborOffset;
		}
	    }
    

	    // TODO: implement a "goto-neighbor" operation
	    // yielding a vertex_iterator! -> useful for 
	    // watershed algo.


	    CoordsGridGraphNeighborIterator & operator++()
	    {
		++index_;
		const shape_type& neighborOffset((*neighborhood_)[index_]);
		tmp_ = pos_ + neighborOffset;
		return *this;
	    }

	    CoordsGridGraphNeighborIterator  operator++(int)
	    {
		CoordsGridGraphNeighborIterator ret(*this);
		++*this;
		return ret;
	    }

	    reference // FIXME: reference to temporary to avoid copies?
	    operator*()
	    {
		return tmp_;
	    }

	    const_reference 
	    operator*() const
	    {
		return tmp_;
	    }
    
	    size_t index() const
	    {
		return index_;
	    }

	    size_t neighborIndex() const
	    {
		return (*neighborIndexLookup_)[index_];
	    }
    
	    bool operator==(self_type const & other) const
	    {
		return index_ == other.index();
	    }
    
	    bool operator!=(self_type const & other) const
	    {
		return index_ != other.index();
	    }

	    bool hasMore() {
		return index_ != (*neighborIndexOffsets_).size();
	    }

	protected:
	    // FIXME: This could possibly be contracted to a single 
	    //   vertex_descriptor element,
	    //   with op++ using as something like pos_ += offset[newindex] - offset[oldindex] 
	    vertex_descriptor pos_;
	    vertex_descriptor tmp_;
	    size_t index_;
	    const ArrayVector<shape_type> *neighborhood_;
	    const ArrayVector<MultiArrayIndex> *neighborIndexOffsets_;
	    const ArrayVector<int> *neighborIndexLookup_;
	};



    } // namespace detail
} // namespace vigra



namespace std {
    template<unsigned int N>
    ostream& operator<<(ostream& out,
			const vigra::detail::CoordsGridGraphNeighborIterator<N>& arg)
    {
	out << "nb" << arg.index();
	return out;
    }
};

#endif // VIGRA_MULTI_GRIDGRAPH_COORDS_NEIGHBOR_ITERATOR_HXX
