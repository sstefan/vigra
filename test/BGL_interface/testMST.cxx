// 3D tests:
#define TEST_3D_REFERENCE_IMPL
//#define TEST_3D_REFERENCE_IMPL2 // vigra-side part broken
#define TEST_LISTBASED
#define BGL_COMPAT // requires TEST_LISTBASED
#define TEST_BGL_COMPAT  // requires BGL_COMPAT
#define TEST_BGL_GRIDGRAPH_COORDS  // requires BGL_COMPAT
#define TEST_BGL_GRIDGRAPH_PTRS  // requires BGL_COMPAT


#include <typeinfo>
#include <iostream>
#include <list>

#include <unittest.hxx>
#include <vigra/imageinfo.hxx> // this includes multi_iterator.hxx 
#include <vigra/impex.hxx>
#include <vigra/copyimage.hxx>
#include <vigra/voxelneighborhood.hxx>
#include <vigra/localminmax.hxx>

#include <vigra/timing.hxx>
#include <vigra/functorexpression.hxx>

#ifdef TEST_LISTBASED
// need to include BGL-related bits these before the algorithms:
#include <vigra/multi_gridgraph_coords.hxx>
#endif


// for multi-dim image experiments:
#include <vigra/multi_array.hxx>
#include <vigra/multi_localminmax.hxx>

#ifdef TEST_BGL_COMPAT
#ifdef WITH_BOOST_GRAPH
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#endif
#include <vigra/graph_algorithms.hxx>
#endif



//#define TEST_SIMILARITY_KERNEL exp(Param(-1.0/4)*sq((Arg1()-Arg2())))
#define TEST_SIMILARITY_KERNEL abs(Arg1()-Arg2())



USETICTOC

using namespace std;
using namespace vigra;

const bool verbose = true;
const size_t repetitions = 1;
const size_t best_of_repetitions = 3;


typedef vigra::MultiArray<2, unsigned char> image_2D_type;
typedef MultiArrayShape<2>::type Shp2D;


struct TestImageReader {

    typedef vigra::MultiArray<2, unsigned char> image_2D_type;
    typedef MultiArrayShape<2>::type Shp2D;
    typedef vigra::MultiArrayView<2, image_2D_type::value_type> multiarrayview_type;

    image_2D_type in;
    image_2D_type out;

    multiarrayview_type v_src;
    multiarrayview_type v_dest;

    TestImageReader() 
    {
	const string input = "test.png";

	// read image using vigra
	if (!isImage(input.c_str())) {
	    throw runtime_error("Test input image not readable.");
	}
	cout << "Reading image " << input << endl;
	vigra::ImageImportInfo info(input.c_str());
	
	Shp2D imsize(info.size());
	cout << "INPUT IMAGE SIZE: " << imsize << endl;

	in = image_2D_type(imsize);
	out = image_2D_type(imsize);
	
	if (!info.isGrayscale()) {
	    vigra::BRGBImage intmp(info.width(), info.height());
	    importImage(info, destImage(intmp));
	    // convert to grayscale (luminance channel):
	    copyImage(srcIterRange(intmp.upperLeft(), intmp.lowerRight(), 
				   RGBToGrayAccessor<vigra::BRGBImage::PixelType>()), 
		      destImage(in));
	    
	} else {
	    importImage(info, destImage(in));
	}

	// get a MultiArrayView on the input image
	multiarrayview_type v_src = multiarrayview_type(in);
	multiarrayview_type v_dest = multiarrayview_type(out);
	cout << "SOURCE VIEW SIZE:" << v_src.shape() << endl;
    }

    void writeOutput(string output="out.png") { // FIXME: add some basename string
	// string output("out.png");
	cout << "Writing result " << output << endl;
	exportImage(srcImageRange(out), vigra::ImageExportInfo(output.c_str()));
    }
};


struct TestData_2D 
{
    enum { dim=2 };
    typedef vigra::MultiArray<dim, double> multiarray_type;
    typedef multiarray_type::value_type value_type;

    multiarray_type vol_src, vol_dest;
    TestData_2D() {
	multiarray_type::size_type sz;
	sz[0] = 5;
	sz[1] = 2;
	vol_src = multiarray_type(sz);
	double data[] = {
	    8,7,1,4,9,
	    9,5,0,2,5
	};
	const double *i = data;
	for (multiarray_type::iterator iter=vol_src.begin(); iter != vol_src.end(); ++iter, ++i) {
	    *iter = *i;
	}
	vol_src = vol_src.transpose();
    }

    template<class GRAPH, class PREDECESSORS>
    bool
    checkPredecessors(const GRAPH &graph, 
		      const PREDECESSORS& predecessors, 
		      bool indirectNeighbors,
		      bool verbose=false) 
    {
	bool res = true;
	using namespace vigragraph;
	typedef typename vigragraph::graph_traits<GRAPH>::vertex_descriptor v;

	if (indirectNeighbors) {
	    res &= get(predecessors, v(0,0)) == v(0,0);  // start node
	    res &= get(predecessors, v(1,0)) == v(0,0);
	    res &= get(predecessors, v(0,1)) == v(0,0);
	    res &= get(predecessors, v(1,1)) == v(0,1);
	    res &= get(predecessors, v(0,2)) == v(1,1);
	    res &= get(predecessors, v(1,2)) == v(0,2);
	    res &= get(predecessors, v(0,3)) == v(1,3);
	    res &= get(predecessors, v(1,3)) == v(0,2);
	    res &= get(predecessors, v(0,4)) == v(1,4);
	    res &= get(predecessors, v(1,4)) == v(0,3);
	} else {
	    res &= get(predecessors, v(0,0)) == v(0,0);  // start node
	    res &= get(predecessors, v(1,0)) == v(0,0);
	    res &= get(predecessors, v(0,1)) == v(0,0);
	    res &= get(predecessors, v(1,1)) == v(0,1);
	    res &= get(predecessors, v(0,2)) == v(1,2);
	    res &= get(predecessors, v(1,2)) == v(1,1);
	    res &= get(predecessors, v(0,3)) == v(1,3);
	    res &= get(predecessors, v(1,3)) == v(1,2);
	    res &= get(predecessors, v(0,4)) == v(1,4);
	    res &= get(predecessors, v(1,4)) == v(1,3);
	}

	return res;
    }

    template<class GRAPH, class EDGELIST>
    bool
    checkMSTEdgeList(const GRAPH &graph, 
		     const EDGELIST& edgelist, 
		      bool indirectNeighbors,
		      bool verbose=false) 
    {
	bool res = true;
	using namespace vigragraph;
	typedef typename vigragraph::graph_traits<GRAPH>::vertex_descriptor v;

	// check (using AdjacencyMatrix access)
	// if all required MST edges are in the list (and also not more)

	if (indirectNeighbors) {
	    res &= edgelist.size() == 9;
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,0), v(1,0), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,0), v(0,1), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,1), v(1,1), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(1,1), v(0,2), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,2), v(1,2), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,2), v(1,3), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,3), v(1,3), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,3), v(1,4), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,4), v(1,4), graph).first);
	} else {
 	    res &= edgelist.size() == 9;
 	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,0), v(1,0), graph).first);
 	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,0), v(0,1), graph).first);
 	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,1), v(1,1), graph).first);


	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(1,2), v(1,1), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(1,1), v(1,2), graph).first); // same edge

	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,2), v(1,2), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(1,2), v(1,3), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,3), v(1,3), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(1,3), v(1,4), graph).first);
	    res &= edgelist.end() != std::find(edgelist.begin(), edgelist.end(), edge(v(0,4), v(1,4), graph).first);
	}

	return res;
    }
};


// {{{ TestData_3D

// 3D test data made from stacking a 2D input image.
// FIXME:  add suitable check
struct TestData_3D 
{
    enum { dim=3 };
    typedef vigra::MultiArray<3, image_2D_type::value_type> multiarray_type;
    typedef multiarray_type::view_type multiarrayview3d_type;
    //   typedef vigra::MultiArrayView<3, vigra::BImage::PixelType> multiarrayview3d_type;
    typedef vigra::MultiArrayView<2, image_2D_type::value_type> multiarrayview_type;

    typedef multiarray_type::value_type value_type;

    TestImageReader testImageReader;

    const image_2D_type &in_;
    image_2D_type out;
    multiarray_type vol_src, vol_dest;
    multiarray_type::size_type sz;


    /// create a 3D test graph out of a 2D image
    // (alternatingly stack the test image with zero-images into a test volume)
    TestData_3D()
    // const image_2D_type &in)
	: testImageReader(), in_(testImageReader.in), out(in_.shape())
    {
	multiarrayview_type v_src(in_);

	sz[2]=20;
	sz[0]=in_.shape()[0];
	sz[1]=in_.shape()[1];

	vol_src = multiarray_type(sz, 0.0);
	for (size_t z=0; z<20; z+=2) {
	    vol_src.bind<2>(z).copy(v_src);
	}
    }
  
    bool
    check(bool verbose=false) 
    {
	bool res = true;
	// FIXME provide result check
	return res;
    }

    template<class GRAPH, class PREDECESSORS>
    bool
    checkPredecessors(const GRAPH &graph, 
		      const PREDECESSORS& prec, 
		      bool indirectNeighbors,
		      bool verbose=false) 
    {
	return true; // cannot check this currently
    }
};



template<class T>
void save_volume(const vigra::MultiArrayView<3, T> &view, std::string fnbase = std::string("slice")) 
{
    for (size_t z=0; z<view.shape()[2]; ++z) {
	std::stringstream fn;
	fn << fnbase << z << ".png";
	//	exportImage(srcImageRange(view.bind<2>(z)), vigra::ImageExportInfo(fn.str().c_str())); // clang doesn't understand this
 	exportImage(srcImageRange(view.bindOuter(z)), vigra::ImageExportInfo(fn.str().c_str()));
    }
}


// }}}


namespace vigragraph {
    namespace helpers {
#ifdef TEST_BGL_COMPAT
	template<class GRAPH>
	struct DummyEdgeWeightMap : vigragraph::put_get_helper<const double&, DummyEdgeWeightMap<GRAPH> >
	{
	    typedef edge_property_tag type;

	    typedef double value_type;
	    typedef const double& reference;
	    typedef readable_property_map_tag category;

	    typedef typename graph_traits<GRAPH>::edge_descriptor key_type;

	    value_type constantvalue;
	    DummyEdgeWeightMap(value_type constantvalue_ = 1.0) : constantvalue(constantvalue_) {}

	    inline
	    reference operator[](const key_type &index) const {
		return constantvalue;
	    }
	};


	// An ad-hoc computed edge_weight property map computed from the node weights
	// Take a kernel (functor) as parameter.
	// TODO: Functor not necessarily constant?!
	template<class GRAPH, class VERTEXWEIGHTMAP, class FUNCTOR>
	struct PixelSimilarityEdgeWeightMap : get_helper<const double, PixelSimilarityEdgeWeightMap<GRAPH,VERTEXWEIGHTMAP,FUNCTOR> >
	{
	    typedef edge_property_tag type;
	    typedef readable_property_map_tag category;

	    typedef typename graph_traits<GRAPH>::edge_descriptor key_type;
	    typedef double value_type;
	    typedef const double& reference;

	    const GRAPH &graph_;
	    const VERTEXWEIGHTMAP &inputmap_;
	    const FUNCTOR &functor_;

	    PixelSimilarityEdgeWeightMap(const GRAPH &graph,
					 const VERTEXWEIGHTMAP &inputmap,
					 const FUNCTOR &functor) : graph_(graph), inputmap_(inputmap), functor_(functor) {
	    }

	    inline
	    value_type operator[](const key_type &index) const {
		// using namespace vigragraph;
		value_type tmp = functor_(get(inputmap_, target(index, graph_)), get(inputmap_, source(index, graph_)));
// 		cout << "PM: for input "<<get(inputmap_, target(index, graph_))<<","<<get(inputmap_, source(index, graph_))
// 		     <<  " at " << target(index, graph_) << " and " << source(index, graph_) 
// 		     << " computed value " << tmp << endl;
		return tmp;
	    }
	};

	// corresponding simple creator function:
	template<class GRAPH, class VERTEXWEIGHTMAP, class FUNCTOR>
	PixelSimilarityEdgeWeightMap<GRAPH, VERTEXWEIGHTMAP, FUNCTOR>
	make_PixelSimilarityEdgeWeightMap(const GRAPH &g,
					  const VERTEXWEIGHTMAP &inputmap,
					  const FUNCTOR &functor) {
	    return PixelSimilarityEdgeWeightMap<GRAPH, VERTEXWEIGHTMAP, FUNCTOR>(g, inputmap, functor);
	}

	
	template<class GRAPH>
	static size_t countEdges(const GRAPH& graph, bool verbose=false) {
	    typedef typename graph_traits<GRAPH>::edge_iterator edge_iterator;
	    edge_iterator i,ie;
	    tie(i, ie) = edges(graph);
	    size_t count = 0;
	    for (; i != ie; ++i) {
		++count;
		if (verbose) {
		    std::cout << "Edge " << count
			      << ": edge_descriptor=" << *i 
			      << " (from " << source(*i, graph) << " to " << target(*i, graph) << ")"
			      << std::endl;
		}
	    }
	    return count;
	}

#endif
    } // namespace helpers
} // namespace vigragraph

struct MSTTests {
#ifdef TEST_BGL_COMPAT


#if 1
    template<class TESTDATA>
    void testND_BGL_Prim_CoordsVariant(bool directNeighborhood=false)
    { 
	//typedef TestData_3D TestData;
	typedef TESTDATA TestData;
	TestData td;
	const unsigned int dim = TestData::dim;
	cout << endl << "Testing BGL Prim's MST algorithm on vigra grid graph; Coords variant, " << (directNeighborhood? "":"in")<<"direct neighborhood variant." << endl;
	using namespace vigragraph;

	typedef typename TestData::value_type value_type;
	typedef GridGraphView_CoordsDescriptor<dim> grid_view_type;
	typedef typename graph_traits<grid_view_type>::vertex_descriptor vertex_descriptor;
	grid_view_type ggv(td.vol_src.shape(), directNeighborhood ? DirectNeighborhood : IndirectNeighborhood);

	using namespace vigra::functor;

	typedef MultiArrayView_property_map<typename TestData::multiarray_type::view_type> grid_data_pm_type;
	grid_data_pm_type weights(td.vol_src);

#if 1
	if (1) {
	    typedef typename graph_traits<grid_view_type>::edge_descriptor edge_descriptor;

	    // typedef helpers::DummyEdgeWeightMap<grid_view_type> edge_weight_map_type;
	    //	typedef PixelSimilarityEdgeWeightMap<grid_view_type, grid_data_pm_type> edge_weight_map_type;
	    // FIXME: there must be some nicer way to denote the type of a functor expression like this one:
	    //       -> there is, using a creator function; however, we then cannot have this typedef...

	    typedef helpers::PixelSimilarityEdgeWeightMap<grid_view_type, grid_data_pm_type, typeof(TEST_SIMILARITY_KERNEL)> edge_weight_map_type;
	    edge_weight_map_type edgeWeights(ggv, weights,  TEST_SIMILARITY_KERNEL);

	    // an explicitly stored edge_weight map
	    // TODO: need some convenience method to construct such an edge property map
	    typedef MultiArray<dim+1, double> edge_weight_storage_type;
	    typedef MultiArrayView_property_map<typename edge_weight_storage_type::view_type> edge_weight_propmap_type;
	    edge_weight_storage_type explicitEdgeStorage(ggv.edge_propmap_shape());
	    edge_weight_propmap_type explicitEdgeWeights(explicitEdgeStorage);
	    cout << " EDGE PM SHAPE:" << explicitEdgeStorage.shape() << endl;
	    helpers::copy_edge_property_map(ggv, edgeWeights, explicitEdgeWeights);
	    // compare_propmaps(graph, edgeWeights, explicitEdgeWeights);
	    // ... run MST on explicit variant ... compare results as well ... TODO
	}
#endif	

	// required for algorithm output:
	// predecessor map: use an image of size_t's here...	
	//	typedef	MultiArray<dim, MultiArrayIndex> predecessor_array_type;
	typedef	MultiArray<dim, vertex_descriptor> predecessor_array_type;
	predecessor_array_type predecessor_data(td.vol_src.shape(), vertex_descriptor(0.0));
	MultiArrayView_property_map<typename predecessor_array_type::view_type> predecessors(predecessor_data);


	// test copying vertex property maps:
	helpers::copy_vertex_property_map(ggv, predecessors, predecessors);


#if WITH_BOOST_GRAPH
	// access markers/weights "property maps" via corresponding get/set functions
	{
	    cout << "Timing Prim MST Test " << dim << "D (new, list-based, via BGL (vertex_desc=coord), direct propmaps src/dest): " << endl;
	    cout << "GRID GRAPH SIZE: " << prod(ggv.shape()) << endl;
	    //	    TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)

	    // required concepts for prim's algorithm:
	    // Vertex List Graph and Incidence Graph. -> easier as that is readily available 

	    //	    boost::prim_minimum_spanning_tree(ggv, predecessors, boost::weight_map(edgeWeights));
	    boost::prim_minimum_spanning_tree(ggv, predecessors, 
					      boost::weight_map(helpers::make_PixelSimilarityEdgeWeightMap(ggv, weights, TEST_SIMILARITY_KERNEL)));

	    // output if small enough:
	    if (prod(ggv.shape()) <= 20) {
		// iterate vertices, output predecessors:
		typedef typename vigragraph::graph_traits<grid_view_type>::vertex_iterator vertex_iterator;
		vertex_iterator i,ie;
		vigragraph::tie(i, ie) = vigragraph::vertices(ggv);
		for (; i != ie; ++i) {
		    cout << get(weights, *i) << " in " << *i << " <- " << get(predecessors, *i) << endl;
		}
	    }
	    // to visualize result: read out predecessor map and denote the edge orientation in a result image...
	    // or just plot the resulting graph in 3D (e.g. using VRML? X3DOM? WebGL?)
	    should(td.checkPredecessors(ggv, predecessors, !directNeighborhood, true));
	}
#endif // WITH_BOOST_GRAPH
    }
#endif    


    template<class TESTDATA>
    void testND_BGL_Prim_CoordsVariant_IndirectNeighborhood()
    {
	testND_BGL_Prim_CoordsVariant<TESTDATA>(false);
    }

    template<class TESTDATA>
    void testND_BGL_Prim_CoordsVariant_DirectNeighborhood()
    {
	testND_BGL_Prim_CoordsVariant<TESTDATA>(true);
    }


#if 1
    template<class TESTDATA>
    void testND_BGL_Prim_CoordsVariant_DirectNeighborhoodCleaned()
    { 
	typedef TESTDATA TestData;
	TestData td;
	const unsigned int dim = TestData::dim;
	cout << "Testing BGL Prim's MST algorithm on vigra grid graph; Coords variant, direct neighborhood." << endl;
	using namespace vigragraph;

	typedef typename TESTDATA::value_type value_type;
	typedef GridGraphView_CoordsDescriptor<dim> grid_view_type;
	typedef typename grid_view_type::vertex_descriptor vertex_descriptor;
	grid_view_type ggv(td.vol_src.shape(), true /*=directNeighborhood*/);

	using namespace vigra::functor;

	typedef MultiArrayView_property_map<typename TestData::multiarray_type::view_type> grid_data_pm_type;
	grid_data_pm_type weights(td.vol_src);
	// required for algorithm output:
	// predecessor map: use an image of size_t's here...	
	//	typedef	MultiArray<dim, MultiArrayIndex> predecessor_array_type;
	typedef	MultiArray<dim, vertex_descriptor> predecessor_array_type;
	predecessor_array_type predecessor_data(td.vol_src.shape(), vertex_descriptor(0.0));
	MultiArrayView_property_map<typename predecessor_array_type::view_type> predecessors(predecessor_data);


#if WITH_BOOST_GRAPH
	// access markers/weights "property maps" via corresponding get/set functions
	{
	    boost::prim_minimum_spanning_tree(ggv, predecessors, 
					      boost::weight_map(helpers::make_PixelSimilarityEdgeWeightMap(ggv, weights, TEST_SIMILARITY_KERNEL)));

	    // output if small enough:
	    if (prod(ggv.shape()) <= 20) {
		// iterate vertices, output predecessors:
		typedef typename vigragraph::graph_traits<grid_view_type>::vertex_iterator vertex_iterator;
		vertex_iterator i,ie;
		vigragraph::tie(i, ie) = vigragraph::vertices(ggv);
		for (; i != ie; ++i) {
		    cout << get(weights, *i) << " in " << *i << " <- " << get(predecessors, *i) << endl;
		}
	    }
	    should(td.checkPredecessors(ggv, predecessors, false, true));
	}
#endif // WITH_BOOST_GRAPH
    }
#endif    



    template<class TESTDATA>
    void testND_BGL_Kruskal_CoordsVariant(bool directNeighborhood = false)
    { 
#if 1
	cout << "Testing BGL Kruskal's MST algorithm on vigra grid graph; Coords variant, " << (directNeighborhood? "":"in")<<"direct neighborhood variant." << endl;
	typedef TESTDATA TestData;
	TestData td;
	using namespace vigragraph;
	const unsigned int dim = TestData::dim;

	typedef typename TestData::value_type value_type;
	typedef GridGraphView_CoordsDescriptor<dim> grid_view_type;
	grid_view_type ggv(td.vol_src.shape(), directNeighborhood ? DirectNeighborhood : IndirectNeighborhood);

	typedef MultiArrayView_property_map<typename TestData::multiarray_type::view_type> grid_data_pm_type;
	grid_data_pm_type weights(td.vol_src);

	typedef typename graph_traits<grid_view_type>::edge_descriptor edge_descriptor;


#if WITH_BOOST_GRAPH
	using namespace vigra::functor;

	typedef helpers::PixelSimilarityEdgeWeightMap<grid_view_type, grid_data_pm_type, typeof(TEST_SIMILARITY_KERNEL)> edge_weight_map_type;
	edge_weight_map_type edgeWeights(ggv, weights,  TEST_SIMILARITY_KERNEL);

	// access markers/weights "property maps" via corresponding get/set functions
	{
	    cout << "Timing Kruskal Test 3D (new, list-based, via BGL (vertex_desc=coord), direct propmaps src/dest): " << endl;
	    //	    TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)

	    // required concepts for prim's algorithm:
	    // Vertex List Graph and Incidence Graph. -> easier as that is readily available 
	    

	    // CHECK: Why necessary to pass the vertex_index map directly?
	    typename property_map<grid_view_type, vertex_index_t>::type
		vertexIndex = get(vigragraph::vertex_index, ggv);

	    // a Kruskal example:
	    std::list<edge_descriptor> mst;
	    boost::kruskal_minimum_spanning_tree(ggv, 
						 std::back_inserter(mst), 
						 weight_map(edgeWeights).
						 vertex_index_map(vertexIndex));
  
	    for(typename std::list<edge_descriptor>::iterator it = mst.begin(); it != mst.end(); ++it)
		{
		    std::cout << " edge in MST: " << *it 
			      << " from " << source(*it, ggv)
			      << " to " << target(*it, ggv)
			      << std::endl;
		}
	    should(td.checkMSTEdgeList(ggv, mst, !directNeighborhood, true));
	}
#endif // WITH_BOOST_GRAPH
#endif    
    }


    template<class TESTDATA>
    void testND_BGL_Kruskal_CoordsVariant_IndirectNeighborhood()
    {
	testND_BGL_Kruskal_CoordsVariant<TESTDATA>(false);
    }

    template<class TESTDATA>
    void testND_BGL_Kruskal_CoordsVariant_DirectNeighborhood()
    {
	testND_BGL_Kruskal_CoordsVariant<TESTDATA>(true);
    }



    // another test case that popped up:
    // check the edge descriptor "reversal" for undirected edge descriptor uniqueness
    void testEdgeReversal() 
    {
	using namespace vigragraph;
	typedef TestData_2D TestData;
	bool directNeighborhood = true;
	TestData td;
	const unsigned int dim = TestData::dim;
	typedef GridGraphView_CoordsDescriptor<dim> Graph;
	Graph graph(td.vol_src.shape(), directNeighborhood ? DirectNeighborhood : IndirectNeighborhood);
	typedef graph_traits<Graph>::vertex_descriptor v;
	
	shouldEqual(edge(v(1,1), v(1,2), graph).first, edge(v(1,2), v(1,1), graph).first);
	shouldEqual(edge(v(1,1), v(1,2), graph).second, edge(v(1,2), v(1,1), graph).second);
    }


    // test different ways to iterate graph's edges
    // and index into an edge property map
    void testEdgePropertyMaps() 
    {
	using namespace vigragraph;
	typedef TestData_2D TestData;
	bool directNeighborhood = true;
	TestData td;
	const unsigned int dim = TestData::dim;
	typedef GridGraphView_CoordsDescriptor<dim> Graph;
	Graph graph(td.vol_src.shape(), directNeighborhood ? DirectNeighborhood : IndirectNeighborhood);

	property_map<Graph, vertex_index_t>::type
	    vertexIndex = get(vertex_index, graph);

	// An explicitly stored undirected edge_weight map:
	// a standard MultiArray is used as property map via adapter

	// TODO: need some convenience method to construct such an edge property map
	typedef MultiArray<dim+1, pair<size_t, size_t> > edge_weight_storage_type;
	typedef MultiArrayView_undirected_edge_property_map<edge_weight_storage_type::view_type, Graph> edge_weight_propmap_type;
	edge_weight_storage_type explicitEdgeStorage(graph.edge_propmap_shape());
	explicitEdgeStorage.init(pair<size_t, size_t>(0,0));
	edge_weight_propmap_type explicitEdgeWeights(explicitEdgeStorage, graph);
	cout << " UNDIRECTED EDGE PM SHAPE:" << explicitEdgeStorage.shape() << endl;

	// iterate vertices, then out_edges
	typedef vigragraph::graph_traits<Graph>::vertex_iterator vertex_iterator;
	typedef vigragraph::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
	vertex_iterator v,ve;
	out_edge_iterator oe,oee;
	vigragraph::tie(v, ve) = vigragraph::vertices(graph);
	size_t count = 0;
	for (; v != ve; ++v) {
	    vigragraph::tie(oe, oee) = vigragraph::out_edges(*v, graph);
	    for (; oe != oee; ++oe) {
		// store source and target vertex_descriptor SOI in edge map
		++count;
		put(explicitEdgeWeights, *oe, pair<size_t, size_t>(get(vertexIndex, source(*oe, graph)),
								   get(vertexIndex, target(*oe, graph))));
	    }
	}
	
	// iterate via edge_iterator and check stored values
	typedef vigragraph::graph_traits<Graph>::edge_iterator edge_iterator;
	edge_iterator ei,eie;
	size_t count2 = 0;
	vigragraph::tie(ei, eie) = vigragraph::edges(graph);
	for (; ei != eie; ++ei) {
	    ++count2;
	    should((get(explicitEdgeWeights, *ei) ==
		    pair<size_t, size_t>(get(vertexIndex, source(*ei, graph)),
					 get(vertexIndex, target(*ei, graph))))
		   || 
		   (get(explicitEdgeWeights, *ei) ==
		    pair<size_t, size_t>(get(vertexIndex, target(*ei, graph)),
					 get(vertexIndex, source(*ei, graph)))			   
		    ));
	}
	shouldEqual(count, 2*count2);
	// TODO: possibly also check via edge(v1, v2) and vice versa via source, target
    }


    void testND_BGL_MaxFlow_CoordsVariant()
    { 
#if 0
	cout << "Testing BGL Boykov Kolmogorov Max-Flow on vigra grid graph; Coords variant." << endl;
	// requirements for graph:
	// * directed (yet to be implemented!) VertexListGraph,EdgeListGraph,IncidenceGraph.
	//   (all the reverse edges need to be in the graph)

	// Edmonds-Karp needs only a (directed) VertexListGraph / IncidenceGraph.

	// ... we would anyways need to construct some wrapper graph for the grid graph that includes the s,t nodes.

	// TODO
#endif    
    }


    void testAdjacencyMatrixAccess()
    {
	using namespace vigragraph;
	// 2D grids 
	if (1) {
	    typedef GridGraphView_CoordsDescriptor<2> Graph;
	    typedef Graph::vertex_descriptor v;
	    typedef Graph::edge_descriptor e;
	    
	    Graph::shape_type size(4, 3);
	    Graph graphDirect(size, DirectNeighborhood);
	    Graph graphIndirect(size, IndirectNeighborhood);

	    // check edge existence:
	    shouldEqual(edge(v(0,0),v(0,1),graphDirect).second, true);
	    shouldEqual(edge(v(0,1),v(0,0),graphDirect).second, true);
	    shouldEqual(edge(v(0,0),v(0,1),graphIndirect).second, true);
	    shouldEqual(edge(v(0,1),v(0,0),graphIndirect).second, true);

	    shouldEqual(edge(v(0,0),v(0,0),graphDirect).second, false);
	    shouldEqual(edge(v(0,0),v(0,0),graphIndirect).second, false);
	    shouldEqual(edge(v(1,1),v(1,1),graphDirect).second, false);
	    shouldEqual(edge(v(1,1),v(1,1),graphIndirect).second, false);

	    shouldEqual(edge(v(0,0),v(1,1),graphDirect).second, false);
	    shouldEqual(edge(v(1,1),v(0,0),graphDirect).second, false);
	    shouldEqual(edge(v(0,0),v(1,1),graphIndirect).second, true);
	    shouldEqual(edge(v(1,1),v(0,0),graphIndirect).second, true);

	    shouldEqual(edge(v(0,0),v(2,0),graphDirect).second, false);
	    shouldEqual(edge(v(0,0),v(2,0),graphIndirect).second, false);
	    shouldEqual(edge(v(0,0),v(2,2),graphIndirect).second, false);

	    // check edge endpoints
	    {
		e e1 = edge(v(0,0),v(1,0), graphDirect).first;
		cout << "Failed Edge: " << e1 << endl;
		should((source(e1, graphDirect) == v(0,0)) || (target(e1, graphDirect) == v(0,0)));
		should((source(e1, graphDirect) == v(1,0)) || (target(e1, graphDirect) == v(1,0)));
		should(source(e1, graphDirect) != target(e1, graphDirect));
	    }
	    {
		e e1 = edge(v(0,0),v(1,0), graphDirect).first;
		should((source(e1, graphDirect) == v(0,0)) || (target(e1, graphDirect) == v(0,0)));
		should((source(e1, graphDirect) == v(1,0)) || (target(e1, graphDirect) == v(1,0)));
		should(source(e1, graphDirect) != target(e1, graphDirect));
	    }
	    {
		e e1 = edge(v(0,0),v(1,0), graphIndirect).first;
		should((source(e1, graphIndirect) == v(0,0)) || (target(e1, graphIndirect) == v(0,0)));
		should((source(e1, graphIndirect) == v(1,0)) || (target(e1, graphIndirect) == v(1,0)));
		should(source(e1, graphIndirect) != target(e1, graphIndirect));
	    }
	    {
		e e1 = edge(v(0,1),v(0,0), graphDirect).first;
		should((source(e1, graphDirect) == v(0,0)) || (target(e1, graphDirect) == v(0,0)));
		should((source(e1, graphDirect) == v(0,1)) || (target(e1, graphDirect) == v(0,1)));
		should(source(e1, graphDirect) != target(e1, graphDirect));
	    }
	    {
		e e1 = edge(v(0,1),v(0,0), graphIndirect).first;
		should((source(e1, graphIndirect) == v(0,0)) || (target(e1, graphIndirect) == v(0,0)));
		should((source(e1, graphIndirect) == v(0,1)) || (target(e1, graphIndirect) == v(0,1)));
		should(source(e1, graphIndirect) != target(e1, graphIndirect));
	    }
	    {
		e e1 = edge(v(0,0),v(1,1), graphIndirect).first;
		should((source(e1, graphIndirect) == v(0,0)) || (target(e1, graphIndirect) == v(0,0)));
		should((source(e1, graphIndirect) == v(1,1)) || (target(e1, graphIndirect) == v(1,1)));
		should(source(e1, graphIndirect) != target(e1, graphIndirect));
	    }
	    {
		e e1 = edge(v(1,1),v(0,0), graphIndirect).first;
		should((source(e1, graphIndirect) == v(0,0)) || (target(e1, graphIndirect) == v(0,0)));
		should((source(e1, graphIndirect) == v(1,1)) || (target(e1, graphIndirect) == v(1,1)));
		should(source(e1, graphIndirect) != target(e1, graphIndirect));
	    }
	}
    }


    void testNumEdges() 
    {
	// The direct neighborhood edge count can be checked with Mathematica's
	// EdgeCount function, e.g. EdgeCount[GridGraph[{3, 2, 2, 2}]] yields 52

	using namespace vigragraph;

	// 1D grid
	if (1) {
	    typedef GridGraphView_CoordsDescriptor<1> Graph;
	    Graph::shape_type size(5);
	    Graph graphDirect(size, DirectNeighborhood);
	    Graph graphIndirect(size, IndirectNeighborhood);

 	    shouldEqual(num_edges(graphDirect), 4);
	    shouldEqual(num_edges(graphIndirect), 4);
	    // check by counting:
	    shouldEqual(helpers::countEdges(graphDirect), 4);
	    shouldEqual(helpers::countEdges(graphIndirect), 4);
	}


	// 2D grids
	if (1) {
	    typedef GridGraphView_CoordsDescriptor<2> Graph;
	    
	    Graph::shape_type size(2,2);
	    Graph graphDirect(size, DirectNeighborhood);
	    Graph graphIndirect(size, IndirectNeighborhood);

 	    shouldEqual(num_edges(graphDirect), 4);
 	    shouldEqual(num_edges(graphIndirect), 6);
	    // check by counting:
	    shouldEqual(helpers::countEdges(graphDirect), 4);
	    shouldEqual(helpers::countEdges(graphIndirect), 6);
	}


	// 2D grids 
	if (1) {
	    typedef GridGraphView_CoordsDescriptor<2> Graph;
	    
	    Graph::shape_type size(4, 3);
	    Graph graphDirect(size, DirectNeighborhood);
	    Graph graphIndirect(size, IndirectNeighborhood);
	    cout << "GRAPH SHAPE:" << graphDirect.shape() << endl;
	    shouldEqual(helpers::countEdges(graphDirect), 9+8); // 2xy-2x-2y
	    shouldEqual(helpers::countEdges(graphIndirect), 9+8+6+6);
 	    shouldEqual(num_edges(graphDirect), 9+8);
 	    shouldEqual(num_edges(graphIndirect), 9+8+6+6);
	}
	// 3D grids
	if (1) {
	    typedef GridGraphView_CoordsDescriptor<3> Graph;
	    int x=3, y=2, z=2;
	    Graph::shape_type size(x, y, z);
	    Graph graphDirect(size, DirectNeighborhood);
	    Graph graphIndirect(size, IndirectNeighborhood);

	    shouldEqual(helpers::countEdges(graphDirect), 3*x*y*z - x*y - x*z - y*z); 
	    shouldEqual(helpers::countEdges(graphIndirect), 
			(x-1)*y*z + x*(y-1)*z + x*y*(z-1) + // 3 axis-aligned edge orientations
			2*(x-1)*(y-1)*z + 2*x*(y-1)*(z-1) + 2*(x-1)*y*(z-1) + // 6 face-diagonal edge orientations
			4*(x-1)*(y-1)*(z-1) // 4 space diagonal orientations
			);
 	    shouldEqual(num_edges(graphDirect), 20);
 	    shouldEqual(num_edges(graphIndirect), 50);
	}

	// 4D grids
	if (1) {
	    typedef GridGraphView_CoordsDescriptor<4> Graph;
	    int x=3, y=2, z=2, t=2;
	    Graph::shape_type size(x, y, z, t);
	    Graph graphDirect(size, DirectNeighborhood);
	    Graph graphIndirect(size, IndirectNeighborhood);

	    shouldEqual(helpers::countEdges(graphDirect), 52);
	    shouldEqual(helpers::countEdges(graphIndirect), 212);  // closed formula still missing...
 	    shouldEqual(num_edges(graphDirect), 52);
            shouldEqual(num_edges(graphIndirect), 212);
	}


	// 5D grids
	if (1) {
	    typedef GridGraphView_CoordsDescriptor<5> Graph;
	    int x=3, y=2, z=2, t=2, u=4;
	    Graph::shape_type size(x, y, z, t, u);
	    Graph graphDirect(size, DirectNeighborhood);
	    Graph graphIndirect(size, IndirectNeighborhood);

	    shouldEqual(helpers::countEdges(graphDirect), 280);
	    shouldEqual(helpers::countEdges(graphIndirect), 2192);  // closed formula still missing...
 	    shouldEqual(num_edges(graphDirect), 280);
            shouldEqual(num_edges(graphIndirect), 2192);
	}

    }



    void checkConcepts_CoordsVariant() {
#if WITH_BOOST_GRAPH
	cout << "Testing BGL concepts on grid graph; Coords variant." << endl;
	using namespace vigragraph;
	typedef GridGraphView_CoordsDescriptor<3> grid_graph_type;

	BOOST_CONCEPT_ASSERT((GraphConcept<grid_graph_type>));
	BOOST_CONCEPT_ASSERT((VertexListGraphConcept<grid_graph_type>));
	BOOST_CONCEPT_ASSERT((IncidenceGraphConcept<grid_graph_type>));
	BOOST_CONCEPT_ASSERT((EdgeListGraphConcept<grid_graph_type>));
	BOOST_CONCEPT_ASSERT((AdjacencyMatrixConcept<grid_graph_type>));

	// TODO:
	// Requirements for the temporary property maps for the algorithm:
	//  rank_map, predecessor_map; or alternatively provide a vertex_index_map to use the defaults.

// 	typedef MultiArrayView_property_map<TestData_3D::multiarrayview3d_type> data_propmap_type;
// 	BOOST_CONCEPT_ASSERT((ReadablePropertyMapConcept<data_propmap_type, grid_graph_type::vertex_descriptor>));
// 	BOOST_CONCEPT_ASSERT((WritablePropertyMapConcept<data_propmap_type, grid_graph_type::vertex_descriptor>));
// 	BOOST_CONCEPT_ASSERT((LvaluePropertyMapConcept<data_propmap_type, grid_graph_type::vertex_descriptor>));
#endif
    }

#endif // TEST_BGL_COMPAT


};




struct GridGraphsTestSuite
    : public vigra::test_suite
{
    GridGraphsTestSuite()
	: vigra::test_suite("GridGraphsTestSuite")
    {

#ifdef TEST_BGL_COMPAT
	// concept checks:
	add( testCase( &MSTTests::checkConcepts_CoordsVariant));
#ifdef TEST_BGL_GRIDGRAPH_COORDS
 	add( testCase( &MSTTests::testEdgeReversal));
 	add( testCase( &MSTTests::testEdgePropertyMaps));
 	add( testCase(( &MSTTests::testAdjacencyMatrixAccess )));
 	add( testCase(( &MSTTests::testNumEdges )));
 	add( testCase(( &MSTTests::testND_BGL_Prim_CoordsVariant_IndirectNeighborhood<TestData_2D>))); 
  	add( testCase(( &MSTTests::testND_BGL_Prim_CoordsVariant_DirectNeighborhood<TestData_2D>))); 
//  	//add( testCase(( &MSTTests::testND_BGL_Prim_CoordsVariant<TestData_3D>))); 
 	add( testCase(( &MSTTests::testND_BGL_Kruskal_CoordsVariant_IndirectNeighborhood<TestData_2D>))); 
 	add( testCase(( &MSTTests::testND_BGL_Kruskal_CoordsVariant_DirectNeighborhood<TestData_2D>))); 
#endif

#endif

    }
};


int main(int argc, char ** argv) {
    GridGraphsTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);

    return 0;
}



/// Local Variables: 
/// c-basic-offset: 4
/// End:
