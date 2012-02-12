#ifndef VIGRA_MULTI_GRIDGRAPH_COORDS_OUT_EDGE_ITERATOR_HXX
#define VIGRA_MULTI_GRIDGRAPH_COORDS_OUT_EDGE_ITERATOR_HXX

#include "graphs.hxx"
#include "tinyvector.hxx"

namespace vigra {
    namespace detail {

	template<unsigned int N>
	class CoordsGridGraphOutEdgeIterator
	{
	public:
	    typedef typename MultiArrayShape<N>::type shape_type;

	    typedef shape_type vertex_descriptor;
	    typedef typename MultiArrayShape<N+1>::type edge_descriptor;

	    typedef typename CoupledScanOrderIteratorFactory<N>::coupled_iter_type vertex_iterator;

	    typedef CoordsGridGraphOutEdgeIterator<N> self_type;

	    typedef edge_descriptor value_type;

	    typedef value_type*  pointer;
	    typedef value_type const * const_pointer;
	    typedef value_type& reference;
	    typedef value_type const & const_reference;

	    typedef std::ptrdiff_t difference_type;
	    typedef std::forward_iterator_tag  iterator_category;


	    CoordsGridGraphOutEdgeIterator() : index_(0),
					       neighborIndexLookup_(0)
	    {}

	    CoordsGridGraphOutEdgeIterator(vertex_descriptor pos,
					   const ArrayVector<shape_type> &neighborhood,
					   const ArrayVector<MultiArrayIndex> &neighborIndexOffsets,
					   const ArrayVector<int> &neighborIndexLookup,
					   bool end=false)
		: pos_(pos),
		  index_(0),
		  neighborIndexLookup_(&neighborIndexLookup)
	    {
		if (end)
		    index_ += neighborIndexLookup.size();
		else {
		    TinyVectorView<typename edge_descriptor::value_type, N>(tmp_.data()) = pos;
		    tmp_[N] = neighborIndex();
		}
	    }
    

	    CoordsGridGraphOutEdgeIterator & operator++()
	    {
		++index_;
		tmp_[N] = neighborIndex();
		return *this;
	    }

	    CoordsGridGraphOutEdgeIterator  operator++(int)
	    {
		CoordsGridGraphOutEdgeIterator ret(*this);
		++*this;
		return ret;
	    }

	    reference 
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
		return index_ != (*neighborIndexLookup_).size();
	    }

	protected:
	    vertex_descriptor pos_; // Maybe eliminate later
	    edge_descriptor tmp_;
	    size_t index_;
	    const ArrayVector<int> *neighborIndexLookup_;
	};


    } // namespace detail
} // namespace vigra


#endif // VIGRA_MULTI_GRIDGRAPH_COORDS_OUT_EDGE_ITERATOR_HXX
