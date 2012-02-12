#ifndef VIGRA_MULTI_GRIDGRAPH_COORDS_EDGE_ITERATOR_HXX
#define VIGRA_MULTI_GRIDGRAPH_COORDS_EDGE_ITERATOR_HXX


#include "graphs.hxx"
#include "tinyvector.hxx"

namespace vigra {
    namespace detail {


	// Edge iterator for undirected graphs. 
	// Composed of a vertex_iterator and an out_edge_iterator
	// (which in case of the undirected graph is filtered for the unique edges).


	template<class GRAPH>
	class CoordsGridGraphEdgeIterator
	{
	public:
	    typedef CoordsGridGraphEdgeIterator<GRAPH> self_type;
	    typedef typename vigragraph::graph_traits<GRAPH>::vertex_iterator vertex_iterator;
	    typedef typename vigragraph::graph_traits<GRAPH>::vertex_descriptor vertex_descriptor;
	    typedef typename vigragraph::graph_traits<GRAPH>::out_edge_iterator out_edge_iterator;
	    typedef typename vigragraph::graph_traits<GRAPH>::edge_descriptor edge_descriptor;

	    typedef  edge_descriptor value_type;

	    typedef value_type*  pointer;
	    typedef value_type const * const_pointer;
	    typedef value_type& reference;
	    typedef value_type const & const_reference;

	    typedef std::ptrdiff_t difference_type;
	    //    typedef std::input_iterator_tag  iterator_category;
	    typedef std::forward_iterator_tag  iterator_category;


	    CoordsGridGraphEdgeIterator() : graph_(0)
	    {}

	    CoordsGridGraphEdgeIterator(const GRAPH &graph);


	    CoordsGridGraphEdgeIterator & inc();
    
	    CoordsGridGraphEdgeIterator & operator++();

	    CoordsGridGraphEdgeIterator  operator++(int)
	    {
		CoordsGridGraphEdgeIterator ret(*this);
		++*this;
		return ret;
	    }

	    reference 
	    operator*()
	    {
		return *outEdgeIterator_;
	    }

	    const_reference 
	    operator*() const
	    {
		return *outEdgeIterator_;
	    }
    
	    bool operator==(self_type const & other) const
	    {
		if (graph_!=other.graph_)
		    return false;
		if ((graph_ != 0) && (other.graph_ != 0)) {
		    if (vertexIterator_ != other.vertexIterator_) 
			return false;
		    if (outEdgeIterator_ != other.outEdgeIterator_) 
			return false;
		}
		return true;
	    }
    
	    bool operator!=(self_type const & other) const
	    {
		return !operator==(other);
	    }

	protected:
	    vertex_iterator vertexIterator_,vertexIteratorEnd_,;
	    out_edge_iterator outEdgeIterator_,outEdgeIteratorEnd_;
	    const GRAPH *graph_;
	};




    } // namespace detail
} // namespace vigra


#endif // VIGRA_MULTI_GRIDGRAPH_COORDS_EDGE_ITERATOR_HXX
