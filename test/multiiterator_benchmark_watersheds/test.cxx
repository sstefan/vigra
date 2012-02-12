#define TEST_WATERSHEDS_3D_REF
// #define TEST_WATERSHEDS_3D
// #define TEST_WATERSHEDS_4D

#define BGL_COMPAT // requires TEST_LISTBASED
#define TEST_BGL_COMPAT  // requires BGL_COMPAT
#define TEST_BGL_GRIDGRAPH_COORDS  // requires BGL_COMPAT
// #define TEST_BGL_GRIDGRAPH_PTRS  // requires BGL_COMPAT


 #define WRITE_TEST_VOLUMES // write out data for debugging; if undefined, save_volume is a dummy call
#define SIDELENGTH_3D 50
#define SIDELENGTH_4D 20

#ifdef WITH_BOOST_GRAPH
#define NONGRID
#endif

#include <typeinfo>
#include <iostream>

#include <unittest.hxx>

#include <vigra/imageinfo.hxx> // this includes multi_iterator.hxx 
#include <vigra/impex.hxx>

#include <vigra/timing.hxx>

// for multi-dim image experiments:
#include <vigra/multi_array.hxx>


#ifdef TEST_BGL_COMPAT  // requires BGL_COMPAT
#ifdef TEST_BGL_GRIDGRAPH_COORDS  // requires BGL_COMPAT
#include <vigra/multi_gridgraph_coords.hxx>
#endif
#ifdef TEST_BGL_GRIDGRAPH_PTRS  // requires BGL_COMPAT
#include <vigra/multi_gridgraph_ptroffset.hxx>
#endif
#endif


#include "multi_testdata.hxx"


#ifdef NONGRID
#include "graphtest.hxx"
#endif


// Watersheds:
#ifdef TEST_WATERSHEDS_3D_REF
#include <vigra/multi_pointoperators.hxx>
#include <vigra/watersheds3d.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/seededregiongrowing3d.hxx>
#endif



USETICTOC

using namespace std;
using namespace vigra;

const bool verbose = true;
const size_t repetitions = 1;

#ifdef WRITE_TEST_VOLUMES
const size_t best_of_repetitions = 1;
#else
const size_t best_of_repetitions = 3;
#endif

struct TestImageReader {

    typedef vigra::MultiArray<2, unsigned char> image_2D_type;
    typedef MultiArrayShape<2>::type Shp2D;

    image_2D_type in;
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
    }
};



struct WatershedTests : TestImageReader {
};


struct WatershedTestsSyntheticImage {

    //    typedef unsigned char  ValueType; // test gives nice "ringing" effect on larger images when overflowing
    typedef double  ValueType;
    typedef unsigned char  QuantizedValueType;

    void watershedTestReferenceImplementation() {
#ifdef TEST_WATERSHEDS_3D_REF
	// Watershed Tests:
	// (adapted from Vigra example for 2D)
	{
	    cout << "Testing Watershed code:" << endl;
	    typedef testdata::TestData3_nD<ValueType, 3> TestData;
	    TestData td(SIDELENGTH_3D);
	    testdata::save_volume(td.vol_src, "testvol");
	    // a labeling image of same size
	    typedef vigra::MultiArray<3, unsigned int> labelarray_type;
	    labelarray_type labeling(td.sz, 0.0);
        
	    // compute a gradient magnitude volume for input:
	    // TODO; for now, copy image itself
	    // TestData_3D::multiarrayview3d_type gradMag(td.vol_src);
	    TestData::view_type gradMag(td.vol_src);
	    unsigned int max_region_label;

	    cout << "Timing 3D SRG(reference impl.): " << endl;
	    TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
		// compute seeds beforehand (use connected components of all pixels 
		// where the gradient  is below 4.0)
		max_region_label = 
		generateWatershedSeeds3D(srcMultiArrayRange(gradMag), 
					 destMultiArray(labeling),
					 NeighborCode3DTwentySix(),
					 //NeighborCode3DSix(),
					 //WatershedOptions().seed_options.levelSets(4.0));
	                                 SeedOptions().levelSets(4.0));
	    cout << "Max region label: " << max_region_label << endl;
	    testdata::save_volume(labeling, "labeling_init");

	    // quantize the gradient image to 256 gray levels
	    typedef vigra::MultiArray<3, QuantizedValueType> quantized_type;
	    quantized_type gradMag256(td.vol_src.shape(), static_cast<QuantizedValueType>(0));

	    vigra::FindMinMax<float> minmax; 
	    inspectMultiArray(srcMultiArrayRange(gradMag), minmax); // find original range
	    transformMultiArray(srcMultiArrayRange(gradMag), destMultiArray(gradMag256),
				linearRangeMapping(minmax, 0, 255));
	    testdata::save_volume(gradMag256, "gradMag256");

	    // create a statistics functor for region growing
	    detail::WatershedStatistics<quantized_type::value_type, unsigned int> regionstats;

	    seededRegionGrowing3D(srcMultiArrayRange(gradMag256), 
				  srcMultiArray(labeling),
				  destMultiArray(labeling), 
				  regionstats, 
				  KeepContours,
				  NeighborCode3DTwentySix());
	    TICTOCLOOP_END
#ifdef WRITE_TEST_VOLUMES
		cout << "  (Including file I/O)" << endl;
#endif
	    // cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
        
	    // FIXME: The 2D variant of watershedsRegionGrowing includes the seed computation???
	    // (only when option "mini" is not set to "Unspecified".)

	    // 	  // call the turbo algorithm with 256 bins, using 8-neighborhood
	    // 	  watershedsRegionGrowing(srcImageRange(gradMag256), destImage(labeling),
	    // 				  WatershedOptions().turboAlgorithm(256));
	    cout << "Max region label after region growing: " << max_region_label << endl;
	    testdata::save_volume(labeling, "labeling_res");

	}
#endif
    }

    void watershedTest3D() 
    {
#ifdef TEST_WATERSHEDS_3D
	{
	    const int dim=3;
	    typedef TestData3_nD<ValueType, dim> TestData;
	    cout << "Testing Watershed code:" << endl;
	    // FIXME: Use a different (synthetic!) test image here, 
	    // which actually has local 3D minima!
	    TestData td(SIDELENGTH_3D);
	    //	testdata::save_volume(td.vol_src, "testvol");
	    // a labeling image of same size
	    typedef vigra::MultiArray<dim, unsigned int> labelarray_type;
	    labelarray_type labeling(td.sz, 0.0);
        
	    // compute a gradient magnitude volume for input:
	    // TODO; for now, copy image itself
	    TestData::view_type gradMag(td.vol_src);
	    unsigned int max_region_label;

	    cout << "watershedsRegionGrowing_nD():" << endl;
	    max_region_label = 
		watershedsRegionGrowing_nD(gradMag,
					   labeling,
					   NeighborCode3DTwentySix(), // FIXME generalize!
					   WatershedOptions()); // FIXME: Test different options
	    cout << "Max region label after watershedsRegionGrowing_nD: " << max_region_label << endl;
	}
#endif
    }


    void watershedTest4D()  
    {
#ifdef TEST_WATERSHEDS_4D
      {
	  const int dim=4;
	  typedef TestData3_nD<ValueType, dim> TestData;
	  cout << "Testing Watershed code:" << endl;
	  TestData td(SIDELENGTH_4D);
	  // a labeling image of same size
	  typedef vigra::MultiArray<dim, unsigned int> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  TestData::view_type gradMag(td.vol_src);
	  unsigned int max_region_label;

	  cout << "watershedsRegionGrowing_nD():" << endl;
	  max_region_label = 
	      watershedsRegionGrowing_nD(gradMag,
					 labeling,
 					 NeighborCode3DTwentySix(), // FIXME generalize!
					 WatershedOptions()); // FIXME: Test different options
	  cout << "Max region label after watershedsRegionGrowing_nD: " << max_region_label << endl;
      }
#endif
    }


    void watershedTest3D_BGL_coords()  
    {
#if defined(TEST_BGL_COMPAT) && defined(TEST_BGL_GRIDGRAPH_COORDS)
      {
	  using namespace vigragraph;
	  const int dim=3;
	  typedef ValueType ValueType;
	  typedef unsigned int LabelType;

	  typedef testdata::TestData3_nD<ValueType, dim> TestData;
	  cout << "Testing Watershed code via BGL gridgraph(coords variant):" << endl;
	  TestData td(SIDELENGTH_3D);
	  // a labeling image of same size
	  typedef vigra::MultiArray<dim, LabelType> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  TestData::view_type gradMag(td.vol_src);
	  unsigned int max_region_label;


	  typedef GridGraphView_CoordsDescriptor<dim> grid_graph_type;
	  // grid_graph_type graph(gradMag.shape(), DirectNeighborhood);
	  grid_graph_type graph(gradMag.shape(), IndirectNeighborhood);

	  // wrap the data as external property maps compatible with this graph:
	  MultiArrayView_property_map<TestData::view_type> gradMagMap(gradMag);
	  MultiArrayView_property_map<labelarray_type::view_type> destMap(labeling);


	  // algorithm below requires a temporary regions array
	  // (check if it can be avoided!)
	  typedef vigra::MultiArray<dim, int> tmpRegions_type;
	  tmpRegions_type tmpRegions(td.sz, 0.0);
	  MultiArrayView_property_map<tmpRegions_type::view_type> tmpRegionsMap(tmpRegions);


	  // for now, test seededregiongrowing directly instead of using the 
	  // watershed wrapper:
	  // destMap is used as seed as well.

	  const double max_cost = 0.0;
	  // const double max_cost = 1e100;
	  // create a statistics functor for region growing
	  vigra::detail::WatershedStatistics<ValueType, LabelType> regionstats;

	  cout << "watershedsRegionGrowing_gridgraph():" << endl;
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions);
	  // firstly, generate seeds:
#if 0
	  max_region_label = 
	    generateWatershedSeeds_graph(graph,
					 gradMagMap, 
					 destMap, 
					 tmpRegionsMap, // used as temporary labeling array
					 WatershedOptions().seed_options);
#else
	  // threshold variant for seed generation:
	  max_region_label = 
	    generateWatershedSeeds_graph(graph,
					 gradMagMap, 
					 destMap, 
					 tmpRegionsMap, // used as temporary labeling array
					 WatershedOptions().seed_options.levelSets(4.0));
	  // write destMap:
	  cout << "Max region label after generateWatershedSeeds_gridgraph: " << max_region_label << endl;
	  testdata::save_volume(labeling, "gg_labeling_init");
#endif

	  // second, run region growing:
	  seededRegionGrowing_graph(graph,
				    gradMagMap,
				    destMap,
				    destMap,
				    tmpRegionsMap,
				    regionstats,
				    // SRGType(CompleteGrow),
				    SRGType(KeepContours),
				    max_cost
				    );
	  TICTOCLOOP_END
	  cout << "Max region label after watershedsRegionGrowing_gridgraph: " << max_region_label << endl;
	  testdata::save_volume(labeling, "gg_labeling_res");
      }
#endif
    }



    void watershedTest3Dwrapper_BGL_coords()  
    {
#if defined(TEST_BGL_COMPAT) && defined(TEST_BGL_GRIDGRAPH_COORDS)
      {
	  using namespace vigragraph;
	  const int dim=3;
	  typedef ValueType ValueType;
	  typedef unsigned int LabelType;

	  typedef testdata::TestData3_nD<ValueType, dim> TestData;
	  cout << "Testing Watershed code via BGL gridgraph(coords variant):" << endl;
	  TestData td(SIDELENGTH_3D);
	  // a labeling image of same size
	  typedef vigra::MultiArray<dim, LabelType> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  TestData::view_type gradMag(td.vol_src);
	  unsigned int max_region_label;


	  typedef GridGraphView_CoordsDescriptor<dim> grid_graph_type;
	  // grid_graph_type graph(gradMag.shape(), DirectNeighborhood);
	  grid_graph_type graph(gradMag.shape(), IndirectNeighborhood);

	  // wrap the data as external property maps compatible with this graph:
	  MultiArrayView_property_map<TestData::view_type> gradMagMap(gradMag);
	  MultiArrayView_property_map<labelarray_type::view_type> destMap(labeling);


	  // algorithm below requires a temporary regions array
	  // (check if it can be avoided!)
	  typedef vigra::MultiArray<dim, int> tmpRegions_type;
	  tmpRegions_type tmpRegions(td.sz, 0.0);
	  MultiArrayView_property_map<tmpRegions_type::view_type> tmpRegionsMap(tmpRegions);

	  cout << "direct watershedsRegionGrowing_gridgraph():" << endl;
	  WatershedOptions options;
	  options.keepContours().seed_options.levelSets(4.0);
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions);
	  max_region_label = 
	      watershedsRegionGrowing_graph(graph,
					    gradMagMap,
					    destMap,
					    tmpRegionsMap,
					    options);
	  TICTOCLOOP_END
	  cout << "Max region label after direct call: watershedsRegionGrowing_gridgraph: " << max_region_label << endl;
	  testdata::save_volume(labeling, "gg_labeling_direct_res");
      }
#endif
    }


    void watershedTest4D_BGL_coords()  
    {
#if defined(TEST_BGL_COMPAT) && defined(TEST_BGL_GRIDGRAPH_COORDS)
      {
	  using namespace vigragraph;
	  const int dim=4;
	  typedef ValueType ValueType;
	  typedef unsigned int LabelType;

	  typedef testdata::TestData3_nD<ValueType, dim> TestData;
	  cout << "Testing 4D Watershed code via BGL gridgraph(coords variant):" << endl;
	  TestData td(SIDELENGTH_4D);
	  // a labeling image of same size
	  typedef vigra::MultiArray<dim, LabelType> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  TestData::view_type gradMag(td.vol_src);
	  unsigned int max_region_label;


	  typedef GridGraphView_CoordsDescriptor<dim> grid_graph_type;
	  grid_graph_type graph(gradMag.shape());

	  // wrap the data as external property maps compatible with this graph:
	  MultiArrayView_property_map<TestData::view_type> gradMagMap(gradMag);
	  MultiArrayView_property_map<labelarray_type::view_type> destMap(labeling);


	  // algorithm below requires a temporary regions array
	  // (check if it can be avoided!)
	  typedef vigra::MultiArray<dim, int> tmpRegions_type;
	  tmpRegions_type tmpRegions(td.sz, 0.0);
	  MultiArrayView_property_map<tmpRegions_type::view_type> tmpRegionsMap(tmpRegions);


	  // for now, test seededregiongrowing directly instead of using the 
	  // watershed wrapper:
	  // destMap is used as seed as well.

	  const double max_cost = 0.0;
	  // create a statistics functor for region growing
	  vigra::detail::WatershedStatistics<ValueType, LabelType> regionstats;

	  cout << "watershedsRegionGrowing_gridgraph():" << endl;
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions);
	  // firstly, generate seeds:
	  // local minimum variant not available yet... using threshold instead:
	  max_region_label = 
	    generateWatershedSeeds_graph(graph,
					 gradMagMap, 
					 destMap, 
					 tmpRegionsMap, // used as temporary labeling array
					 WatershedOptions().seed_options.levelSets(4.0));
	  // write destMap:
	  cout << "Max region label after generateWatershedSeeds_gridgraph: " << max_region_label << endl;

	  // second, run region growing:
	  seededRegionGrowing_graph(graph,
				    gradMagMap,
				    destMap,
				    destMap,
				    tmpRegionsMap,
				    regionstats,
				    // SRGType(CompleteGrow),
				    SRGType(KeepContours),
				    max_cost
				    );
	  TICTOCLOOP_END
	  cout << "Max region label after watershedsRegionGrowing_gridgraph: " << max_region_label << endl;
      }
#endif
    }


    void watershedTest4D_BGL_ptrs()  
    {
#if defined(TEST_BGL_COMPAT) && defined(TEST_BGL_GRIDGRAPH_PTRS)
      {
	  const int dim=4;
	  typedef TestData3_nD<ValueType, dim> TestData;
	  cout << "Testing Watershed code via BGL gridgraph(ptrs variant):" << endl;
	  TestData td(SIDELENGTH_4D);
	  // a labeling image of same size
	  typedef vigra::MultiArray<dim, unsigned int> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  TestData::view_type gradMag(td.vol_src);
	  unsigned int max_region_label;

	  cout << "watershedsRegionGrowing_nD():" << endl;
	  max_region_label = 
	      watershedsRegionGrowing_nD(gradMag,
					 labeling,
 					 NeighborCode3DTwentySix(), // FIXME generalize!
					 WatershedOptions()); // FIXME: Test different options
	  cout << "Max region label after watershedsRegionGrowing_nD: " << max_region_label << endl;
      }
#endif
    }


};



#ifdef NONGRID
struct WatershedTestsNonGrid 
    : public testdata::TestGraph {

    void testNonGridWatershed() {
	if (verbose) 
	    cout << "Testing Non-Grid-Graph:" << endl;

	{ 
	    // test a non-grid graph (via BGL)
	    using namespace vigragraph;
	    if (verbose) 
		cout << "Testing Non-Grid-Graph Watershed:" << endl;

	    // create an associated property map taking bools to indicate 

	    // external properties:
	    // TODO: Also test internal property variant for weights!
	    vector<double> &weights_vec(weights_); 
	    vector<int> labels_vec(num_vertices(G));
	    vector<int> tmp_vec(num_vertices(G));

	    // define corresponding property maps:
	    typedef property_map<Graph, vertex_index_t>::type IndexMap;
	    IndexMap index = get(vertex_index, G);

	    // iterator_property_map is responsible for mapping nodes to keys
	    iterator_property_map<vector<double>::iterator, IndexMap>
		weights(weights_vec.begin(), index);
	    iterator_property_map<vector<int>::iterator, IndexMap>
		labels(labels_vec.begin(), index);
	    iterator_property_map<vector<int>::iterator, IndexMap>
		tmp(tmp_vec.begin(), index);

	    // Test with level-set seed computation:
	    {
		WatershedOptions options;
		//	    options.keepContours().seed_options.levelSets(4.0);
		options.keepContours().seed_options.levelSets(2.1);
		int max_region_label = 
		    watershedsRegionGrowing_graph(G,
						  weights,
						  labels,
						  tmp,
						  options);

		cout << "Max region label after watershedsRegionGrowing_gridgraph: " << max_region_label << endl;

		// Output result
		typedef boost::graph_traits<testdata::TestGraph::Graph>::vertex_iterator 
		    graph_scanner;
		graph_scanner it, itend;
		cout << "Graph Watershed Result:" << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << names_[*it];
		}
		cout << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << labels[*it] ;
		}
		cout << endl;
		int expected_labels[] = { 1,1,2,2,2,2,0,2,2,2,2,2,3 };
		int nr=0;
		for (tie(it,itend)=vertices(G); it!=itend; ++it, ++nr) {
		    shouldEqual(expected_labels[nr], labels[*it]);
		}
	    }

	    // Test with local minima seed computation:
	    if (1) {
		setWeights2();
		// write_dot_graph();
		
		WatershedOptions options;
		options.keepContours().seed_options.minima();
		// options.seed_options.minima();
		int max_region_label = 
		    watershedsRegionGrowing_graph(G,
						  weights,
						  labels,
						  tmp,
						  options);
		cout << "Max region label after watershedsRegionGrowing_gridgraph with minima seeds: " << max_region_label << endl;

		// Output result
		typedef boost::graph_traits<testdata::TestGraph::Graph>::vertex_iterator 
		    graph_scanner;
		graph_scanner it, itend;
		cout << "Graph Watershed (seed:localmin) Result:" << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << names_[*it];
		}
		cout << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << labels[*it] ;
		}
		cout << endl;
		// with keepContours():
		// abcdefghijklm
		// 2210222000304
		//without:
		// abcdefghijklm
		// 2212222222324
		int expected_labels[] = { 2, 2, 1, 0, 2, 2, 2, 0, 0, 0, 3, 0, 4 };
		int nr=0;
		for (tie(it,itend)=vertices(G); it!=itend; ++it, ++nr) {
		    shouldEqual(expected_labels[nr], labels[*it]);
		}
	    }


	    // Test with extended local minima seed computation:
	    if (1) {
		setWeights2();
		// TODO: possibly also check seed computation separately:
		//  a b c d e f g h i j k l m
		//  0 0 1 1 0 1 0 0 1 1 1 1 1
		// Seeds after cc-labeling:
		//  a b c d e f g h i j k l m
		//  0 0 1 2 0 3 0 0 2 2 4 2 5

		WatershedOptions options;
		options.keepContours().seed_options.extendedMinima();
		int max_region_label = 
		    watershedsRegionGrowing_graph(G,
						  weights,
						  labels,
						  tmp,
						  options);

		cout << "Max region label after watershedsRegionGrowing_gridgraph with extended minima seeds: " << max_region_label << endl;

		// Output result
		typedef boost::graph_traits<testdata::TestGraph::Graph>::vertex_iterator 
		    graph_scanner;
		graph_scanner it, itend;
		cout << "Graph Watershed (seed:extended localmin) Result:" << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << names_[*it];
		}
		cout << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << labels[*it] ;
		}
		cout << endl;
		int expected_labels[] = {3, 3, 1, 2, 3, 3, 3, 0, 2, 2, 4, 2, 5};
		int nr=0;
		for (tie(it,itend)=vertices(G); it!=itend; ++it, ++nr) {
		    shouldEqual(expected_labels[nr], labels[*it]);
		}
	    }


	    // Test turbo watershed code:
	    if (1) {
		setWeights2();
		// TODO: possibly also check seed computation separately:
		//  a b c d e f g h i j k l m
		//  0 0 1 1 0 1 0 0 1 1 1 1 1
		// Seeds after cc-labeling:
		//  a b c d e f g h i j k l m
		//  0 0 1 2 0 3 0 0 2 2 4 2 5

		WatershedOptions options;
		// FIXME: ADD TEST: using 10 buckets or less should raise exception.
		options.turboAlgorithm(11).seed_options.extendedMinima();
		int max_region_label = 
		    watershedsRegionGrowing_graph(G,
						  weights,
						  labels,
						  tmp,
						  options);

		cout << "Max region label after turbo watershedsRegionGrowing_gridgraph with extended minima seeds: " << max_region_label << endl;

		// Output result
		typedef boost::graph_traits<testdata::TestGraph::Graph>::vertex_iterator 
		    graph_scanner;
		graph_scanner it, itend;
		cout << "Graph Watershed (seed:extended localmin) Result:" << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << names_[*it];
		}
		cout << endl;
		for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		    cout << labels[*it] ;
		}
		cout << endl;
		int expected_labels[] = {3, 3, 1, 2, 3, 3, 3, 1, 2, 2, 4, 2, 5};  // watershed node 'h' now labeled arbitrarily (by 1 here)
		int nr=0;
		for (tie(it,itend)=vertices(G); it!=itend; ++it, ++nr) {
		    shouldEqual(expected_labels[nr], labels[*it]);
		}
	    }

	}
    }

    void testNonGridLocalMax() {
	if (verbose) 
	    cout << "Testing Non-Grid-Graph:" << endl;

	{ 
	    // test a non-grid graph (via BGL)
	    using namespace boost;

	    describe();

	    if (verbose) 
		cout << "Testing Non-Grid-Graph LocalMax:" << endl;

	    // create an associated property map taking bools to indicate 
	    //   if node is locally maximal:

	    // to apply the same localminmax algorithm,
	    // require G to be a model of AdjacencyGraph.
	    // Furthermore, we need a method to associate the boolean label with each node,
	    // i.e. a property map over the nodes of the graph.

	    // external properties:
	    vector<bool> markers_vec(num_vertices(G));
	    vector<double> &weights_vec(weights_); 
	    // TODO: Also test internal property variant for weights!

	    // define corresponding property maps:
	    typedef property_map<Graph, vertex_index_t>::type IndexMap;
	    IndexMap index = get(vertex_index, G);

	    // iterator_property_map is responsible for mapping nodes to keys
	    iterator_property_map<vector<bool>::iterator, IndexMap>
		markers(markers_vec.begin(), index);
	    iterator_property_map<vector<double>::iterator, IndexMap>
		weights(weights_vec.begin(), index);

	    localMinMaxGraph(G,
			     weights,
			     markers,
			     true, // marker
			     (double)0.0, // threshold
			     std::greater<double>());

	    // Output result
	    typedef boost::graph_traits<testdata::TestGraph::Graph>::vertex_iterator 
		graph_scanner;
	    graph_scanner it, itend;
	    cout << "Graph LocalMax Result:" << endl;
	    for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		cout << names_[*it];
	    }
	    cout << endl;
	    for (tie(it,itend)=vertices(G); it!=itend; ++it) {
		cout << (markers[*it] ? "M" : ".");
	    }
	    cout << endl;
	    bool expected_markers[] = { true, false, false, false, true, false, true, 
					false, false, false, false, false, true };
	    int nr=0;
	    for (tie(it,itend)=vertices(G); it!=itend; ++it, ++nr) {
		shouldEqual(expected_markers[nr], markers[*it]);
	    }
	}
    }

};
#endif //  NONGRID



struct GridGraphsTestSuite
: public vigra::test_suite
{
    GridGraphsTestSuite()
    : vigra::test_suite("GridGraphsTestSuite")
    {
	//	add( testCase( &WatershedTests::watershedTestReferenceImplementation));
	add( testCase( &WatershedTestsSyntheticImage::watershedTestReferenceImplementation));
	add( testCase( &WatershedTestsSyntheticImage::watershedTest3D));
	add( testCase( &WatershedTestsSyntheticImage::watershedTest4D));

	add( testCase( &WatershedTestsSyntheticImage::watershedTest3D_BGL_coords)); 
	add( testCase( &WatershedTestsSyntheticImage::watershedTest3Dwrapper_BGL_coords)); 
	add( testCase( &WatershedTestsSyntheticImage::watershedTest4D_BGL_coords)); 
	add( testCase( &WatershedTestsSyntheticImage::watershedTest4D_BGL_ptrs)); 

#ifdef NONGRID
	//	add( testCase( &WatershedTestsNonGrid::testNonGridLocalMax));
	add( testCase( &WatershedTestsNonGrid::testNonGridWatershed));
#endif //  NONGRID

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
