#ifndef GRAPHTEST_HXX
#define GRAPHTEST_HXX

// #include <boost/config.hpp> // why?
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <vector>

#include <boost/graph/graphviz.hpp>

namespace testdata
{

    using namespace std;
    using namespace boost;



    class TestGraph
    {
    public:
	// a directed graph:
	// typedef adjacency_list<listS, vecS> Graph; // VertexList=vecS
	// an undirected graph:
	typedef adjacency_list<listS, vecS, undirectedS> Graph; // VertexList=vecS

	Graph G;

	enum nodenames {a, b, c, d, e, f, g, h, i, j, k, l, m, N};
	vector<string> names_;
	vector<double> weights_;

	typedef iterator_property_map<vector<string>::iterator, property_map<Graph, vertex_index_t>::type>
	  NameMapType;

	typedef iterator_property_map<vector<double>::iterator,  property_map<Graph, vertex_index_t>::type>
	  WeightMapType;

	TestGraph() {
	    string names[] = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"};    
	    names_ = vector<string>(names, names+N);
	    const double weights[] = {2, 1, 2, 6, 5, 3, 5, 4, 6, 3, 3, 5, 1};
	    weights_ = vector<double>(weights, weights+N);

	    G=Graph(N);
	    add_edge(a, b, G);
	    add_edge(b, g, G);
	    add_edge(c, h, G);
	    add_edge(d, i, G);
	    add_edge(e, f, G);
	    add_edge(e, h, G);
	    add_edge(f, g, G);
	    add_edge(g, h, G);
	    add_edge(h, i, G);
	    add_edge(h, k, G);
	    add_edge(i, j, G);
	    add_edge(i, l, G);
	    // m doesn't have any edges
	}

	//! flat weights for testing extendedminima seeds:
	void setWeights2()
	{
	    const double weights[] = {4, 4, 4, 4, 4, 3, 4, 10, 4, 4, 4, 4, 4};
	    weights_ = vector<double>(weights, weights+N);
	};


	NameMapType 
	getNameMap() {
	    return NameMapType(names_.begin(), get(vertex_index, G));
	}

	WeightMapType 
	getWeightMap() {
	    return WeightMapType(weights_.begin(), get(vertex_index, G));
	}

	template <class Name, class Weight>
	class label_writer {
	public:
	    label_writer(Name _name, Weight _weight) : name(_name), weight(_weight) {}
	    template <class VertexOrEdge>
	    void operator()(std::ostream& out, const VertexOrEdge& v) const {
		out << "[label=\"" << name[v] << ":" << weight[v] << "\"]";
	    }
	private:
	    Name name;
	    Weight weight;
	};


	void 
	write_dot_graph() {
	    //	    write_graphviz(std::cout, *this, make_label_writer(getNameMap()));
	    label_writer<NameMapType, WeightMapType>
		lw(getNameMap(),getWeightMap());
	    write_graphviz(std::cout, G, lw);
	}

	void describe() 
	{
	    property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, G);
	    graph_traits<Graph>::vertex_iterator vi, vi_end;
	    graph_traits<Graph>::adjacency_iterator ai, ai_end;
	    for (tie(vi,vi_end) = vertices(G); vi != vi_end; ++vi) {
		cout << names_[get(index_map, *vi)] << " w=" << weights_[get(index_map, *vi)];
		tie(ai, ai_end) = adjacent_vertices(*vi, G);
		if (ai == ai_end)
		    cout << " has no children ";
		else
		    cout << " is parent of ";
		for (; ai != ai_end; ++ai) {
		    cout << names_[get(index_map, *ai)];
		    if (boost::next(ai) != ai_end)
			cout << ", ";
		}
		cout << endl;
	    }

	    // test the out_edge_iterator as well:
	    graph_traits<Graph>::out_edge_iterator oe, oe_end;
	    graph_traits<Graph>::vertex_descriptor 
		s = vertex(0, G);
	    cout << " vertex 0: " << s << endl;

	    for (tie(vi,vi_end) = vertices(G); vi != vi_end; ++vi) {
		cout << "the edges incident to node " <<  names_[get(index_map, *vi)] <<":";

		for (tie(oe, oe_end) = out_edges(*vi, G); oe != oe_end; ++oe)
		    std::cout << " (" << source(*oe, G) 
			      << "," << target(*oe, G) << ")";
		cout << endl;
    
	    }

	} // describe()
    };


}; // Namespace detail

#endif // guard

