namespace vigra {
    namespace detail {

	// Edge iterator for undirected graphs. 
	// Composed of a vertex_iterator and an out_edge_iterator
	// (which in case of the undirected graph is filtered for the unique edges).


	template<class GRAPH>
	CoordsGridGraphEdgeIterator<GRAPH>::CoordsGridGraphEdgeIterator(const GRAPH &graph)
	    : graph_(&graph)
	{
	    vigragraph::tie(vertexIterator_, vertexIteratorEnd_) = vigragraph::vertices(graph);
	    vigragraph::tie(outEdgeIterator_, outEdgeIteratorEnd_) = vigragraph::out_edges(*vertexIterator_, graph);
	    // check if the initial edge needs to be skipped for undirected graph:
	    while (true) {
		if (graph_==0) break;
		if (graph_->map_to_undirected_edge(*outEdgeIterator_) == *outEdgeIterator_) 
		    break;
		inc();
	    }
	}
    

	template<class GRAPH>
	inline
	CoordsGridGraphEdgeIterator<GRAPH> & 
	CoordsGridGraphEdgeIterator<GRAPH>::inc()
	{
	    ++outEdgeIterator_;
	    if (outEdgeIterator_ != outEdgeIteratorEnd_) {
		return *this;
	    }
	    // advance to next node
	    ++vertexIterator_;
	    if (vertexIterator_ != vertexIteratorEnd_) {
		vigragraph::tie(outEdgeIterator_, outEdgeIteratorEnd_) = vigragraph::out_edges(*vertexIterator_, *graph_);
		return *this;
	    }
	    // else set to "invalid" (end) state:
	    graph_ = 0;
	    return *this;
	}
    
	template<class GRAPH>
	inline
	CoordsGridGraphEdgeIterator<GRAPH> & 
	CoordsGridGraphEdgeIterator<GRAPH>::operator++()
	{
	    // filter to skip 'duplicate' edges in undirected graph
	    while (true) {
		inc();
		if (graph_==0) return *this;
		if (graph_->map_to_undirected_edge(*outEdgeIterator_) == *outEdgeIterator_) 
		    break;
		// std::cout << " skipping edge " << *outEdgeIterator_ << std::endl;
	    }
	    return *this;	    
	}

    } // namespace detail
} // namespace vigra
