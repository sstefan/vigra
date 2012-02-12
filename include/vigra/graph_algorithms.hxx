#ifndef VIGRA_GRAPH_ALGORITHMS_HXX
#define VIGRA_GRAPH_ALGORITHMS_HXX

//! some utility functions for graphs

#include <vigra/graphs.hxx>
#include <vigra/functorexpression.hxx>

namespace vigragraph {    

    namespace helpers {
	template<class GRAPH, class COPYFROM, class COPYTO, class FUNCTOR>
	void
	transform_vertex_property_map(const GRAPH &graph, 
				      const COPYFROM &from,
				      COPYTO &to,
				      FUNCTOR functor) {
	    typedef typename graph_traits<GRAPH>::vertex_iterator vertex_iterator;
	    vertex_iterator i,ie;
	    tie(i, ie) = vertices(graph);
	    for (; i != ie; ++i) {
		put(to, *i, functor(get(from, *i)));
	    }
	}

	template<class GRAPH, class COPYFROM, class COPYTO> 
	void
	copy_vertex_property_map(const GRAPH &graph, 
				 const COPYFROM &from,
				 COPYTO &to) {
	    transform_vertex_property_map(graph, from, to, vigra::functor::Arg1());
	}



	template<class GRAPH, class COPYFROM, class COPYTO, class FUNCTOR>
	void
	transform_edge_property_map(const GRAPH &graph, 
				      const COPYFROM &from,
				      COPYTO &to,
				      FUNCTOR functor) {
	    typedef typename graph_traits<GRAPH>::edge_iterator edge_iterator;
	    edge_iterator i,ie;
	    tie(i, ie) = edges(graph);
	    size_t count = 0;
	    for (; i != ie; ++i, ++count) {
		std::cout << " iter edge " << *i << std::endl;
		put(to, *i, functor(get(from, *i)));
	    }
	    std::cout << " PROCESSED " << count << " EDGES." << std::endl;
	}

	template<class GRAPH, class COPYFROM, class COPYTO> 
	void
	copy_edge_property_map(const GRAPH &graph, 
				 const COPYFROM &from,
				 COPYTO &to) {
	    transform_edge_property_map(graph, from, to, vigra::functor::Arg1());
	}

    } // namespace helpers
} // namespace vigragraph

#endif //  VIGRA_GRAPH_ALGORITHMS_HXX
