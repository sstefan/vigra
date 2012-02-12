/************************************************************************/
/*                                                                      */
/*           Copyright 2004-2012 by Ullrich Koethe, Stefan Schmidt      */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */                
/*                                                                      */
/************************************************************************/

//#define VERBOSE


#ifndef VIGRA_MULTI_GRIDGRAPH_COORDS_HXX
#define VIGRA_MULTI_GRIDGRAPH_COORDS_HXX

#include "graphs.hxx"
#include "multi_iterator.hxx"
#include "multi_iterator_coupled.hxx"
#include "mathutil.hxx"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iterator>


#include "multi_gridgraph_neighborhoods.hxx"
#include "multi_gridgraph_coords_neighbor_iterator.hxx"
#include "multi_gridgraph_coords_out_edge_iterator.hxx"
#include "multi_gridgraph_coords_edge_iterator_decl.hxx"


namespace vigra {


//! Grid Graph class to adapt vigra MultiArrayViews to a BGL-like interface.
//       This class only knows about
//       - dimensions
//       - shape
//       - neighborhood type (DirectedNeighborhood or IndirectNeighborhood, i.e.\ including diagonal neighbors)
template<unsigned int N>
class GridGraphView_CoordsDescriptor
{
public:
    typedef GridGraphView_CoordsDescriptor<N> self_type;

    typedef typename MultiArrayShape<N>::type shape_type;
    typedef typename MultiArrayShape<N+1>::type edge_propmap_shape_type;
    typedef MultiArrayIndex index_type;

    typedef typename CoupledScanOrderIteratorFactory<N>::coupled_iter_type vertex_iterator;
    typedef detail::CoordsGridGraphNeighborIterator<N> neighbor_vertex_iterator;
    typedef detail::CoordsGridGraphOutEdgeIterator<N> out_edge_iterator;
    typedef detail::CoordsGridGraphEdgeIterator<self_type> edge_iterator;


    struct traversal_category : virtual public vigragraph::incidence_graph_tag,
                                virtual public vigragraph::adjacency_graph_tag,
                                virtual public vigragraph::vertex_list_graph_tag,
				virtual public vigragraph::edge_list_graph_tag,
                                // virtual public bidirectional_graph_tag,
				virtual public vigragraph::adjacency_matrix_tag
                                { };

    typedef shape_type  vertex_descriptor;
    typedef typename MultiArrayShape<N+1>::type edge_descriptor;

    typedef void in_edge_iterator; // for bidirectional_graph concept, not implemented here

    typedef neighbor_vertex_iterator      adjacency_iterator; // must be a MultiPassInputIterator model

    typedef vigragraph::undirected_tag directed_category;
    typedef vigragraph::disallow_parallel_edge_tag edge_parallel_category;

    typedef MultiArrayIndex     vertices_size_type;
    typedef MultiArrayIndex     edges_size_type;
    typedef MultiArrayIndex     degree_size_type;

    // we only support "external properties".
    typedef vigragraph::no_property vertex_property_type;
    // TODO: Maybe support the vertex -> coordinate map (identity) as the only internal property map
    // and additionally the vertex_descriptor -> ID map (vertex_index = SOI).


    // dummy default cons to satisfy adjacency_graph concept
    GridGraphView_CoordsDescriptor()
    {}


    //! Constructor for grid graph. 
    //  @param shape                  an array of the graph's dimensions as a TinyVector
    //  @param directNeighborsOnly    true for direct neighborhood (axis-aligned edges only) 
    //                                or false for indirect neighborhood (including all diagonal edges)
    GridGraphView_CoordsDescriptor(shape_type const &shape, NeighborhoodType directNeighborsOnly = IndirectNeighborhood) 
	: shape_(shape)
    {
	// use makeArrayNeighborhood to populate neighborhood tables:
	detail::makeArrayNeighborhood(neighborhood, 
				      neighborExists, 
				      causalNeighborhood, 
				      anticausalNeighborhood, 
				      neighborIndexLookup, 
				      directNeighborsOnly);
	
	// compute the neighbor offsets per neighborhood type
	detail::makeArraySubNeighborhood(neighborhood[0], neighborExists, shape_type(1), neighborhoodIndices);

	// compute total number of edges
	num_edges_ = 0;
	for (unsigned int i=0; i<neighborhood[0].size(); ++i) {
	    size_t product = 1;
	    for (unsigned int j=0; j<N; ++j) {
		product *= (neighborhood[0][i][j]==0) ? shape_[j] : shape_[j]-1;
	    }
	    num_edges_ += product;
	}
	num_edges_ /= 2; // because of undirectedness
    }

    inline
    vertex_iterator get_vertex_iterator() const {
	return vertex_iterator(shape_);
    }

    inline
    vertex_iterator get_vertex_end_iterator() const {
	return vertex_iterator(shape_).getEndIterator();
    }

    inline
    neighbor_vertex_iterator get_neighbor_vertex_iterator(const vertex_iterator& pos) const {
	// determine neighborhood type
	unsigned int nbtype = pos.neighborhoodType();
	// instantiate appropriate neighborhood iterator
	return neighbor_vertex_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype]);
    }

    inline
    neighbor_vertex_iterator get_neighbor_vertex_end_iterator(const vertex_iterator &pos) const {
	// determine neighborhood type
	// instantiate appropriate neighborhood iterator end
	unsigned int nbtype = pos.neighborhoodType();
	// instantiate appropriate neighborhood iterator
	return neighbor_vertex_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype], true);
    }


    // --------------------------------------------------
    // support for VertexListGraph:

    inline
    vertices_size_type
    num_vertices() const 
    {
	return prod(shape_);
    }


    // --------------------------------------------------
    // support for IncidenceGraph:

    inline
    out_edge_iterator get_out_edge_iterator(const vertex_iterator& pos) const {
	// determine neighborhood type
	unsigned int nbtype = pos.neighborhoodType();
	// instantiate appropriate neighborhood iterator
	return out_edge_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype]);
    }

    inline
    out_edge_iterator get_out_edge_end_iterator(const vertex_iterator &pos) const {
	// determine neighborhood type
	// instantiate appropriate neighborhood iterator end
	unsigned int nbtype = pos.neighborhoodType();
	// instantiate appropriate neighborhood iterator
	return out_edge_iterator(pos.point(), neighborhood[nbtype], neighborhoodIndices[nbtype], neighborIndexLookup[nbtype], true);
    }

    inline
    const vertex_descriptor& neighborCoordOffset(int fullNeighborhoodIndex) const {
	const shape_type& neighborOffset(neighborhood[0][fullNeighborhoodIndex]);
	return neighborOffset;
    }

    inline 
    degree_size_type out_degree(const vertex_iterator &pos) const {
	// requires to fully reconstructed iterator (to accesss for neighborhood type)
	unsigned int nbtype = pos.neighborhoodType();
	return neighborhood[nbtype].size();
    }


    // --------------------------------------------------
    // support for EdgeListGraph:

    inline
    edges_size_type
    num_edges() const 
    {
	return num_edges_;
    }


    // --------------------------------------------------
    // support for AdjacencyMatrix concept:

    std::pair<edge_descriptor, bool>
    edge(const vertex_descriptor &u, const vertex_descriptor &v) const
    {
	edge_descriptor edge;
	bool found=false;

	vertex_iterator reconstructed = get_vertex_iterator();
	reconstructed += u;
	unsigned int nbtype = reconstructed.neighborhoodType();

	// check if (u-v) in neighborlist (or v-u, to save reconstruction of v!)
	shape_type diff = u-v;
	for (unsigned int i=0; i< neighborhood[nbtype].size(); ++i) {
	    if (diff == neighborhood[nbtype][i]) {
		found = true;
		edge = map_to_undirected_edge(make_edge_descriptor(v, neighborIndexLookup[nbtype][i]));
		break;
	    }
	    else if ((-diff) == neighborhood[nbtype][i]) {
		found = true;
		// need to use edge from v to u in this case:
		edge = map_to_undirected_edge(make_edge_descriptor(u, neighborIndexLookup[nbtype][i]));
		break;
	    }
	}
	return std::make_pair(edge, found);
    }


    // --------------------------------------------------
    // other helper functions:

    inline 
    degree_size_type maxDegree() const {
	// or: return max_degree;
 	return neighborhood[0].size();
    }
    inline 
    degree_size_type halfMaxDegree() const {
 	return maxDegree() / 2;
    }

    inline
    const shape_type& shape() const {
	return shape_;
    }

    static 
    edge_descriptor make_edge_descriptor(const vertex_descriptor &v,
					    index_type nbindex) 
    {
	edge_descriptor res;
	TinyVectorView<typename edge_descriptor::value_type, N>(res.data()) = v;
	res[N] = nbindex;
	return res;
    }

    edge_propmap_shape_type edge_propmap_shape() const {
	edge_propmap_shape_type res;
	TinyVectorView<typename edge_propmap_shape_type::value_type, N>(res.data()) = shape_;
	// res[N] = maxDegree(); // for directed graph
	res[N] = halfMaxDegree(); // for undirected graph
	return res;
    }

    
    //! In case of an undirected graph, the edge u->v is the same as v->u.
    //  This function folds in the edge descriptors corresponding to "causal neighbors",
    //  i.e. those to vertices with lower scan order index, by reversing the edge and computing
    //  the corresponding edge descriptor.
    //  (assuming here the neighbor-indices are in the order of noncausal, causal neighbors....
    //   FIXME: check this again in neighborhood construction!)
    //  (At least, the neighborhood construction is symmetrical, hence this should be OK)
    inline 
    edge_descriptor
    map_to_undirected_edge(const edge_descriptor &e) const {
	edge_descriptor res = e;
	if (res[N] >= halfMaxDegree()) {
	    TinyVectorView<typename edge_descriptor::value_type, N> vertex(res.data());
	    vertex += neighborCoordOffset(res[N]);
	    res[N] = maxDegree() - res[N] - 1;
	}
	return res;
    }



protected:
    shape_type shape_;
    size_t num_edges_;
    ArrayVector<ArrayVector<shape_type> > neighborhood;
    ArrayVector<ArrayVector<MultiArrayIndex> > neighborhoodIndices;
    ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
    ArrayVector<ArrayVector<int> > neighborIndexLookup;
};


    





} // namespace vigra










// Define Traits classes for BGL compatibility:
//   to obtain vertex_iterator, adjacency_iterator etc.

namespace vigragraph
    {
	using namespace vigra;

	template<unsigned int N>
	inline
	std::pair<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator, 
		  typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator >
	vertices(vigra::GridGraphView_CoordsDescriptor<N> &g) 
	{
	    return std::make_pair(g.get_vertex_iterator(),
				  g.get_vertex_end_iterator());    
	}

	// const variant
	template<unsigned int N>
	inline
	std::pair<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator, 
		  typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator >
	vertices(const vigra::GridGraphView_CoordsDescriptor<N> &g) 
	{
	    return std::make_pair(g.get_vertex_iterator(),
				  g.get_vertex_end_iterator());    
	}


	
	template<unsigned int N>
	inline
	typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertices_size_type
	num_vertices(const vigra::GridGraphView_CoordsDescriptor<N> &g) 
	{
	    return g.num_vertices();
	}



	template<unsigned int N>
	inline
	std::pair<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::adjacency_iterator, 
		  typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::adjacency_iterator >
	adjacent_vertices(typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor v,
			  vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{
	    // Here: need to provide a variant that converts the index vertex_descriptor
	    // back into the corresponding node_iterator.
	    // 
	    typedef typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator
		vertex_iterator;
	    vertex_iterator reconstructed = g.get_vertex_iterator();
	    reconstructed += v;
  
	    return std::make_pair(g.get_neighbor_vertex_iterator(reconstructed),
				  g.get_neighbor_vertex_end_iterator(reconstructed));    
	}


	// adjacent_vertices variant in vigra namespace: allows to call adjacent_vertices with vertex_iterator argument
	template<unsigned int N>
	inline
	std::pair<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::adjacency_iterator, 
		  typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::adjacency_iterator >
	adjacent_vertices_at_iterator(typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator const &v,
				      vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{    
	    return std::make_pair(g.get_neighbor_vertex_iterator(v),
				  g.get_neighbor_vertex_end_iterator(v));    
	}



	template<unsigned int N>
	inline
	std::pair<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::out_edge_iterator, 
		  typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::out_edge_iterator >
	out_edges(typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor v,
			  vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{
	    // Here: need to provide a variant that converts the index vertex_descriptor
	    // back into the corresponding node_iterator.
	    // 
	    typedef typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator
		vertex_iterator;
	    vertex_iterator reconstructed = g.get_vertex_iterator();
	    reconstructed += v;
  
	    return std::make_pair(g.get_out_edge_iterator(reconstructed),
				  g.get_out_edge_end_iterator(reconstructed));    
	}


	template<unsigned int N>
	inline
	typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor 
	source(const typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_descriptor e,
			  vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{
	    typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor res =
		TinyVectorView<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_descriptor::value_type, N>(e.data());
	    return res;
	}



	template<unsigned int N>
	inline
	typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor 
	target(const typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_descriptor e,
			  vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{
	    typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor res = 
		TinyVectorView<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_descriptor::value_type, N>(e.data());

	    // the target is a bit more complicated, because we need the help of the graph to find the correct offset:
	    res += g.neighborCoordOffset(e[N]);
	    return res;
	}

	template<unsigned int N>
	inline
	typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::degree_size_type
	out_degree(const typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor v,
			  vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{
	    typedef typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_iterator
		vertex_iterator;
	    vertex_iterator reconstructed = g.get_vertex_iterator();
	    reconstructed += v;
	    return g.out_degree(reconstructed);
	}



	template<unsigned int N>
	inline
	std::pair<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_iterator, 
		  typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_iterator >
	edges(vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{
	    typedef typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_iterator edge_iterator;
	    return std::make_pair(edge_iterator(g), edge_iterator());
	}


	template<unsigned int N>
	inline
	typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edges_size_type
	num_edges(vigra::GridGraphView_CoordsDescriptor<N> const &g) 
	{
	    return g.num_edges();
	}


	// --------------------------------------------------
	// support for AdjacencyMatrix concept:

	template<unsigned int N>
	std::pair<typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::edge_descriptor, bool>
	edge(const typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor &u,
	     const typename vigragraph::graph_traits<vigra::GridGraphView_CoordsDescriptor<N> >::vertex_descriptor &v,
	     vigra::GridGraphView_CoordsDescriptor<N> const &g)
	{
	    return g.edge(u,v);
	}



	// provide get / put for MultiArrayViews, indexed by the above-defined vertex_descriptor (in this case, a coordinate tuple):

	template<class VIEW>
	class MultiArrayView_property_map 
	{
	public:
	    // typedef vigragraph::read_write_property_map_tag category;
	    typedef vigragraph::lvalue_property_map_tag category;
	    //    typedef int value_type;
	    typedef typename VIEW::value_type value_type;
	    typedef typename VIEW::reference reference;
	    typedef typename VIEW::const_reference const_reference;
	    typedef typename vigra::GridGraphView_CoordsDescriptor<VIEW::actual_dimension> graph_type;
	    typedef typename graph_type::vertex_descriptor key_type;
	    MultiArrayView_property_map(const VIEW& view)
		: view_(view) { }
	    template <class T2>
	    inline
	    reference operator[](const T2 & x) { return view_[x]; } 
	    template <class T2>
	    inline
	    const_reference operator[](const T2 & x) const { return view_[x]; } 
    
	protected:
	    VIEW view_;
	};

	template<class VIEW>
	inline
	void put(MultiArrayView_property_map<VIEW> &pmap,
		 const typename MultiArrayView_property_map<VIEW>::key_type &k,
		 const typename MultiArrayView_property_map<VIEW>::value_type& val
		 ) 
	{ 
	    pmap[k] = val;  
	}

	template<class VIEW>
	inline
	typename MultiArrayView_property_map<VIEW>::const_reference 
	get(
	    const MultiArrayView_property_map<VIEW> & pmap, 
	    const typename MultiArrayView_property_map<VIEW>::key_type &k)
	{ 
	    return pmap[k]; 
	}


	//! use a MultiArrayView as an undirected edge property map
	//  (edge_descriptor keys will be transformed by "flipping" the edges if necessary)
	template<class VIEW, class GRAPH>
	class MultiArrayView_undirected_edge_property_map 
	{
	public:
	    // typedef vigragraph::read_write_property_map_tag category;
	    typedef vigragraph::lvalue_property_map_tag category;
	    //    typedef int value_type;
	    typedef typename VIEW::value_type value_type;
	    typedef typename VIEW::reference reference;
	    typedef typename VIEW::const_reference const_reference;
	    typedef GRAPH graph_type;
	    typedef typename graph_type::edge_descriptor key_type;
	    MultiArrayView_undirected_edge_property_map(const VIEW& view, const GRAPH& graph)
		: view_(view), graph_(graph) { }
	    template <class T2>
	    inline
	    reference operator[](const T2 & x) { return view_[graph_.map_to_undirected_edge(x)]; } 
	    template <class T2>
	    inline
	    const_reference operator[](const T2 & x) const { return view_[graph_.map_to_undirected_edge(x)]; } 
    
	protected:
	    VIEW view_;
	    const GRAPH &graph_;
	};

	template<class VIEW, class GRAPH>
	inline
	void put(MultiArrayView_undirected_edge_property_map<VIEW, GRAPH> &pmap,
		 const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::key_type &k,
		 const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::value_type& val
		 ) 
	{ 
	    pmap[k] = val;  
	}

	template<class VIEW, class GRAPH>
	inline
	typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::const_reference 
	get(
	    const MultiArrayView_undirected_edge_property_map<VIEW, GRAPH> & pmap, 
	    const typename MultiArrayView_undirected_edge_property_map<VIEW, GRAPH>::key_type &k)
	{ 
	    return pmap[k]; 
	}




	
	// property map support for mapping coordinates to scan-order indices:

	template<unsigned int N>
	struct IDMapper {
	    typedef typename vigra::GridGraphView_CoordsDescriptor<N> graph_type;
	    typedef vigragraph::readable_property_map_tag category;
	    typedef typename graph_type::index_type value_type;
	    typedef typename graph_type::vertex_descriptor key_type;
	    typedef const value_type& reference;


	    IDMapper(const graph_type &graph) 
		: map_helper(graph.get_vertex_iterator())
	    {}

	    typename graph_type::vertex_iterator map_helper;
	};

	template<unsigned int N>
	struct property_map<vigra::GridGraphView_CoordsDescriptor<N>, vigragraph::vertex_index_t>
	{
	    typedef IDMapper<N> type;
	    typedef IDMapper<N> const_type;
	};


	template<unsigned int N>
	inline
	typename IDMapper<N>::value_type
	get(const IDMapper<N> & mapper, 
	    const typename IDMapper<N>::key_type &k)
	{ 
	    return (mapper.map_helper + k).index();
	}


	template<unsigned int N>
	typename vigragraph::property_map<vigra::GridGraphView_CoordsDescriptor<N>, vigragraph::vertex_index_t>::type
	//typename IDMapper<N>
	get(vigragraph::vertex_index_t, const vigra::GridGraphView_CoordsDescriptor<N> &graph) {
	    // return a lightweight wrapper for the CoupledIterator, which easily allows the conversion of 
	    // coordinates via its += operator followed by index().
	    return IDMapper<N>(graph);
	}

#if 0 
	// CHECK if required: also provide the direct (three-parameter) version for index lookup
	template<unsigned int N>
	typename vigra::GridGraphView_CoordsDescriptor<N>::vertices_size_type
	get(vigragraph::vertex_index_t, 
	    const vigra::GridGraphView_CoordsDescriptor<N> &graph,
	    const typename vigra::GridGraphView_CoordsDescriptor<N>::vertex_descriptor &v) {
	    return (IDMapper<N>(graph).map_helper + v).index();
	}
#endif



	// TODO:
	// eventually provide an edge_index property map as well?
	// (edge_descriptor -> linear contiguous edge index)


} // namespace vigragraph




// FIXME: I'd rather like the definition of this iterator in 
// a single file with the declaration, however I didn't manage
// to resolve the circular dependency between the two templates
// in another way. Have to try some harder sometime else.        
#include "multi_gridgraph_coords_edge_iterator_defn.hxx"



namespace std {
    template<unsigned int N>
    ostream& operator<<(ostream& out,
			const typename vigra::GridGraphView_CoordsDescriptor<N>::vertex_iterator & arg)
    {
	out << "v" << arg.index();
	return out;
    }
};
#endif /* VIGRA_MULTI_GRIDGRAPH_COORDS_HXX */



