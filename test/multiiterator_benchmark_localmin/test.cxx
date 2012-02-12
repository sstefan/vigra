#define TEST_2D
//#define CHECK_NEIGHBORHOODS
#define TEST_2D_ADAPTED_LOCALMINMAX
// 3D tests:
#define TEST_3D_REFERENCE_IMPL
//#define TEST_3D_REFERENCE_IMPL2 // vigra-side part broken
#define TEST_CHAINED_ITERS
#define TEST_LISTBASED
#define TEST_LISTBASED_CHAINED
#define BGL_COMPAT // requires TEST_LISTBASED
#define TEST_BGL_COMPAT  // requires BGL_COMPAT
#define TEST_BGL_INDIRECT_PROPMAPS
#define TEST_BGL_DIRECT_PROPMAPS
#define TEST_BGL_DIRECT_PROPMAPS1
#define TEST_BGL_DIRECT_PROPMAPS2

#include <typeinfo>
#include <iostream>

#include <unittest.hxx>
#include <vigra/imageinfo.hxx> // this includes multi_iterator.hxx 
#include <vigra/impex.hxx>
#include <vigra/copyimage.hxx>
#include <vigra/graphs.hxx>
#include <vigra/voxelneighborhood.hxx>
#include <vigra/localminmax.hxx>

#include <vigra/timing.hxx>

// for multi-dim image experiments:
#include <vigra/multi_array.hxx>
#include <vigra/multi_localminmax.hxx>

#ifdef TEST_LISTBASED
#include <vigra/multi_gridgraph_coords.hxx>
#endif

#ifdef WITH_BOOST_GRAPH
#define TEST_BGL_CONCEPTS
#endif

#ifdef TEST_BGL_CONCEPTS
#include <boost/graph/graph_concepts.hpp>
#endif


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



// {{{ TestData_3D

struct TestData_3D 
{
  typedef vigra::MultiArray<3, image_2D_type::value_type> multiarray_type;
  typedef multiarray_type::view_type multiarrayview3d_type;
//   typedef vigra::MultiArrayView<3, vigra::BImage::PixelType> multiarrayview3d_type;
  typedef vigra::MultiArrayView<2, image_2D_type::value_type> multiarrayview_type;

  typedef multiarray_type::value_type value_type;

  const image_2D_type &in_;
  image_2D_type out;
  const bool allowAtBorder_;
  multiarray_type vol_src, vol_dest;
  multiarray_type::size_type sz;

  /// create a 3D test graph out of a 2D image
  // (alternatingly stack the test image with zero-images into a test volume)
  TestData_3D(const image_2D_type &in, int numNeighbors = 8, bool allowAtBorder = true) 
    : in_(in), out(in.shape()), allowAtBorder_(allowAtBorder)
  {
    // compute reference 2D result:
    vigra::LocalMinmaxOptions options = vigra::LocalMinmaxOptions()
      .neighborhood(numNeighbors)
      .threshold(0)
      .markWith(255);

    if (allowAtBorder)
      options.allowAtBorder();

    vigra::localMaxima(srcImageRange(in), 
		       destImage(out),
		       options);

    multiarrayview_type v_src(in);

    sz[2]=20;
    sz[0]=in.shape()[0];
    sz[1]=in.shape()[1];

    vol_src = multiarray_type(sz, 0.0);
    vol_dest = multiarray_type(sz, 0.0);
    for (size_t z=0; z<20; z+=2) {
      vol_src.bind<2>(z).copy(v_src);
    }
  }
  
  bool
  check(bool verbose=false) 
  {
    const image_2D_type &expected_res_2D(out);
    bool res = true;

    multiarrayview_type v_expected_res(expected_res_2D);
    multiarray_type vol_zeros(sz, 0.0);

    
    // every impair slice should be zero:
    for (size_t z=1; z<20; z+=2) {
      if (vol_dest.bind<2>(z) != vol_zeros.bind<2>(z))
	res = false;
    }
    if (verbose) 
      cout << "impair slices ok?:" << (res?"Yes":"No") << endl;
    // remainder should equal input:
    for (size_t z=0; z<20; z+=2) {
      if (vol_dest.bind<2>(z) != ((z==0) && !allowAtBorder_ ? vol_zeros.bind<2>(z) : v_expected_res)) {
	res = false;
	if (verbose) 
	  cout << "slice " << z << " wrong." << endl;
	// write out slice:
	stringstream fn;
	fn << "slice" << z << ".png";
	exportImage(srcImageRange(vol_dest.bind<2>(z)), vigra::ImageExportInfo(fn.str().c_str()));
      }
    }

    if (!res)
      exportImage(srcImageRange(expected_res_2D), vigra::ImageExportInfo("ref.png"));

    return res;
  }

};

// }}}






struct LocalMinTests : TestImageReader {
    void test2DreferenceImpl() 
    {
      // init output
      out *= 0;
      TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	vigra::localMaxima(
			   srcImageRange(in), 
			   destImage(out),
			   // vigra::LocalMinmaxOptions().neighborhood(4).allowAtBorder().threshold(0).markWith(255)
			   // vigra::LocalMinmaxOptions().neighborhood(8).allowAtBorder().threshold(0).markWith(255)
			   // vigra::LocalMinmaxOptions().neighborhood(8).threshold(0).markWith(255)
			   vigra::LocalMinmaxOptions().neighborhood(4).threshold(0).markWith(255)
			   );
      TICTOCLOOP_END

      writeOutput("out0.png");
    }


#ifdef TEST_2D      
    // these are not really unittests, just checking out the interface...
    void test2D_TraverserTest() 
      {
	cout << "Traverser test: " << endl;
	// tmp: test the multi_iterator: // this iterates one dimension only
	size_t count=0;
	typedef multiarrayview_type::traverser traverser;
	traverser t=v_src.traverser_begin(), tend = v_src.traverser_end();
	for (; t != tend; ++t) {
	  cout << ".";
	  ++count;
	}
	cout << endl << "COUNT=" << count << endl << endl;
      }

    void test2D_StridedScanOrderIteratorTest() 
    {
	cout << "StridedScanOrderIterator test: " << endl;
	// tmp: test the standard scan-order iterator
	size_t count=0;
	typedef multiarrayview_type::iterator scanner;
	scanner t=v_src.begin(), tend = v_src.end();
	for (; t != tend; ++t) {
	  //cout << "T: " << t << endl;
	  //		cout << ",";
	  //	cout << *t;
	  ++count;
	}
	cout << endl << "COUNT=" << count << endl << endl;
    }

#ifdef CHECK_NEIGHBORHOODS
    void test2D_CheckNeighborhoodsTest() 
      // wrap a GridGraphScanOrderIterator around it
      {
	size_t count=0;
	typedef multiarrayview_type::iterator scanner;

	scanner t=v_src.begin(), tend = v_src.end();
	for (; t != tend; ++t) {
	  ++count;

	  scanner::nb_iterator nbit(t.nb1_begin());
	  scanner::nb_iterator nbend(t.nb1_end());
	  
	  if ((count==1) || 
	      (count==2) || 
	      (count==400) || 
	      (count==401) || 
	      (count==402) ||
	      (count==800) ||
	      (count==120000-399) ||
	      (count==120000-398) ||
	      (count==120000)
	      ) {
	    int i=0;
	    for (;nbit != nbend; ++nbit) {
	      ++i;
	      nbit.printDone();
	      cout << " count " << count  << "  NB" <<i << " "<<((void*)&(*nbit)) << endl;
	      //	      break;
	    }
	    cout << endl;
	  }

	  // for testing, output the neighborhood size
	  
	}
	cout << endl << "MyScanner COUNT=" << count << endl << endl;
      }
#endif


#ifdef TEST_2D_ADAPTED_LOCALMINMAX
    void test2D_AdaptedLocalminmaxTest()
    {
      // apply an adapted Version of localminmax
      TIC;
      for (int i=0; i<repetitions; i++)
      {
	// TODO: make function from this loop; generalize, templatify
	double threshold = 0;
	double marker = 255;

	size_t count=0;
	typedef multiarrayview_type::iterator scanner;
	scanner src=v_src.begin(), srcend = v_src.end();
	scanner dest=v_dest.begin(), destend = v_dest.end();
	for (; src != srcend; ++src, ++dest) {
	  const double  refval = *src;
	  if (refval < threshold)
	    continue;
	  ++count;
	  scanner::nb_iterator nbit(src.nb1_begin());
	  scanner::nb_iterator nbend(src.nb1_end());
	  
	  bool local_extremum = true;
	  int j=0;
	  for (;nbit != nbend; ++nbit) {
	    ++j;
	    if (*nbit >= refval) {
	      local_extremum = false;
	      break;
	    }
	  }
	  if (local_extremum)
	    *dest = marker;

	  // TODO: no neighbors => local min/max per definition?
	  
	}
      }
      double time2=TOCN;
      if (verbose) 
	cout << "Timing new: " << repetitions << " repetitions took " << time2 << " ms." << endl;
    }
#endif

    void test2D_LocalMinMaxTest()
    {
      // reinit output
      out *= 0;
      // apply an adapted Version of localminmax
      // TODO
      TIC;
      for (int i=0; i<repetitions; i++)
      {
	localMinMax(v_src, v_dest, multiarrayview_type::value_type(255.0), multiarrayview_type::value_type(0.0), std::greater<multiarrayview_type::value_type>());
      }
      
      double time3=TOCN;
      if (verbose) 
	cout << "Timing new2: " << repetitions << " repetitions took " << time3 << " ms." << endl;

      writeOutput();
    }

#endif // TEST_2D
    

#ifdef CHECK_NEIGHBORHOODS
    void test3D_NeighborhoodCheckTest()
      {
	cout << "3D NEIGHBORHOOD CHECK:" << endl;
	TestData_3D td(in, 4);

	typedef TestData_3D::multiarrayview3d_type view_type;
	typedef view_type::iterator scanner;

	size_t count = 0;
	scanner t=view_type(td.vol_src).begin(), tend = view_type(td.vol_src).end();
	for (; t != tend; ++t) {
	  ++count;
	  scanner::nb_iterator nbit(t.nb1_begin());
	  scanner::nb_iterator nbend(t.nb1_end());
	  
	  if ((count==1) || 
	      (count==2) || 
	      (count==400) || 
	      (count==401) || 
	      (count==402) ||
	      (count==800) ||
	      (count==120000-399) ||
	      (count==120000-398) ||
	      (count==120000)
	      ) {
	    int i=0;
	    for (;nbit != nbend; ++nbit) {
	      ++i;
	      nbit.printDone(); 
	      cout << " count " << count  << "  NB" <<i << " "<<((void*)&(*nbit)) << endl;
	    }
	    cout << endl;
	  }

	}
      }
#endif

#ifdef TEST_3D_REFERENCE_IMPL
    void test3D_ReferenceImplTest()
      {
	cout << "Timing 3D (reference impl localMaxima3D): " << endl;
	TestData_3D td(in, 8, false); // 8-neighborhood in 2D, spare border

        TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	     vigra::localMaxima3D(srcMultiArrayRange(td.vol_src), 
				  destMultiArray(td.vol_dest),
				  255,
				  vigra::NeighborCode3DTwentySix());

	TICTOCLOOP_END
	  cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
      } // 3D test
#endif // TEST_3D_REFERENCE_IMPL

#ifdef TEST_3D_REFERENCE_IMPL2
    // this code has been commented out in Vigra for an unknown reason.
    void test3D_ReferenceImpl2Test()
      {
	cout << "Testing reference impl. (3D, nD multiarray overload version):" << endl;
	TestData_3D td(in, 8, false); // 8-neighborhood in 2D, spare border
	typedef TestData_3D::multiarrayview3d_type view_type;

	TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	       vigra::localMaxima(view_type(td.vol_src), 
				  view_type(td.vol_dest),
				  vigra::LocalMinmaxOptions()
				  .neighborhood(26)
				  .threshold(0)
				  .markWith(255));
	TICTOCLOOP_END
 	 cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
      } // 3D test
#endif // TEST_3D_REFERENCE_IMPL2

#ifdef TEST_CHAINED_ITERS      
    void test3D_ChainedItersTest()
      {
	cout << "Timing 3D (new, chained iterators, 8-neighbors): " << endl;
	//TestData_3D td(in, 8);
	TestData_3D td(in, 4);

	TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	  localMinMax(td.vol_src, td.vol_dest, 
			 TestData_3D::value_type(255.0), 
			 TestData_3D::value_type(0.0), 
			 std::greater<TestData_3D::value_type>());
	TICTOCLOOP_END
 	 cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
      } // 3D test
#endif // TEST_CHAINED_ITERS



#ifdef TEST_LISTBASED

#ifdef TEST_LISTBASED_CHAINED
    void test3D_ListBased_ChainedItersTest()
      { 
	cout << "Timing 3D (new, list-based chained iters): " << endl;
	TestData_3D td(in, 8);

	typedef GridGraphView<3, image_2D_type::value_type> grid_view_type;
	grid_view_type ggv(td.vol_src);

	TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	     localMinMax(ggv, 
			 td.vol_dest, 
			 TestData_3D::value_type(255.0), 
			 TestData_3D::value_type(0.0), 
			 std::greater<TestData_3D::value_type>());
	TICTOCLOOP_END
	 cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
      }
#endif // TEST_LISTBASED_CHAINED


#ifdef TEST_BGL_COMPAT

#ifdef TEST_BGL_INDIRECT_PROPMAPS 
    void test3D_BGL_IndirectPropmapsTest()
      // property map indexing indirect via iterator_property_map
      { 
	cout << "Testing BGL compatibility on grid graph: " << endl;
	TestData_3D td(in, 8);
	using namespace boost;

	typedef image_2D_type::value_type value_type;
	typedef GridGraphView<3, value_type> grid_view_type;
	grid_view_type ggv(td.vol_src);
	grid_view_type ggv_dest(td.vol_dest);

	typedef grid_view_type Graph;

#ifdef TEST_BGL_CONCEPTS
	boost::LvaluePropertyMapConcept<Graph, Graph::vertex_descriptor> checkConcept();

	if (0) { // FIXME: crashes at runtime,
	  //  because a default-constructed graph has no data associated.
	  // Is such concept-checking code supposed to be run at all? 
	  if (verbose) 
	    cout << "concept check:" << endl;
	  // concept checking:
	  AdjacencyGraphConcept<Graph> gcheck;
	  if (verbose) 
	    cout << "concept check OK." << endl;
	}
#endif


	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, ggv);
	//IndexMap index;

	// iterator_property_map is responsible for mapping nodes to keys
	iterator_property_map<grid_view_type::vertex_iterator, IndexMap>
	  markers(ggv_dest.get_vertex_iterator(), index);
	iterator_property_map<grid_view_type::vertex_iterator, IndexMap>
	  weights(ggv.get_vertex_iterator(), index);

	// FIXME: This is what currently makes the code slow according
	//        to the profiler!
	// These iterator_property_maps just need random access to the nodes
	// via a linear index. Is that available cheaper than with the 
	// full-blown StridedScanOrderIterator::operator+=(int)? 
	// (In fact, we do sequential access and not random, thus we don't 
	//  want to pay for the overhead... how is that usually done in BGL?)

	cout << "Timing 3D (new, list-based, via BGL): " << endl;
	TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	localMinMaxGraph(ggv,
			 weights,
			 markers,
			 value_type(255.0), // marker
			 value_type(0.0), // threshold
			 std::greater<value_type>());

	TICTOCLOOP_END
 	 cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
      }
#endif // if TEST_BGL_INDIRECT_PROPMAPS


#ifdef TEST_BGL_DIRECT_PROPMAPS
    void test3D_BGL_DirectPropmapsTest()
      { 
	cout << "Testing BGL compatibility on grid graph; DIRECT ACCESS VARIANT (for source weights)" << endl;
	TestData_3D td(in, 8);
	using namespace boost;

	typedef image_2D_type::value_type value_type;
	typedef GridGraphView<3, value_type> grid_view_type;
	grid_view_type ggv(td.vol_src);
	grid_view_type ggv_dest(td.vol_dest);

	typedef grid_view_type Graph;

	// access markers/weights directly, not via iterator_property_map


	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, ggv);
	//IndexMap index;

	// iterator_property_map is responsible for mapping nodes to keys
	iterator_property_map<grid_view_type::vertex_iterator, IndexMap>
	  markers(ggv_dest.get_vertex_iterator(), index);
//  	iterator_property_map<grid_view_type::vertex_iterator, IndexMap>
//  	  weights(ggv.get_vertex_iterator(), index);

	// just dereference the vertex_iterator (serving as vertex_descriptor)
	// by conversion to its base type? No, that indexes the original graph, 
	// hence it would only work for weights, but not auxiliary ones.
	// Can one do better, by exploiting somehow that the auxiliary property maps
	//  (e.g. markers) usually have the same structure?
	////   scalar_gridgraph_propery_map<...> markers(ggv_dest);

#ifdef TEST_BGL_DIRECT_PROPMAPS1
	{
	  cout << "Timing 3D (new, list-based, via BGL, direct propmaps): " << endl;
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	      localMinMaxGraph(ggv,
			       ggv, // =weights,
			       markers,
			       value_type(255.0), // marker
			       value_type(0.0), // threshold
			       std::greater<value_type>());

	  TICTOCLOOP_END
	  cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
	}
	td.vol_dest *= 0.0; // reset output
#endif // TEST_BGL_DIRECT_PROPMAPS1
#ifdef TEST_BGL_DIRECT_PROPMAPS2
	{
	  cout << "Timing 3D (new, list-based, via BGL, direct propmaps src/dest): " << endl;
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	      localMinMaxGraph2(ggv,
				ggv, // =weights,
				ggv_dest,
				ggv_dest, //=markers,
				value_type(255.0), // marker
				value_type(0.0), // threshold
				std::greater<value_type>());
	  TICTOCLOOP_END
	  cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
	}
#endif // TEST_BGL_DIRECT_PROPMAPS2
      }
#endif // TEST_BGL_DIRECT_PROPMAPS

#endif // TEST_BGL_COMPAT

#endif // TEST_LISTBASED

};




struct GridGraphsTestSuite
: public vigra::test_suite
{
    GridGraphsTestSuite()
    : vigra::test_suite("GridGraphsTestSuite")
    {
	add( testCase( &LocalMinTests::test2DreferenceImpl));
#ifdef TEST_2D      
	add( testCase( &LocalMinTests::test2D_TraverserTest));
	add( testCase( &LocalMinTests::test2D_StridedScanOrderIteratorTest));
  #ifdef CHECK_NEIGHBORHOODS
	add( testCase( &LocalMinTests::test2D_CheckNeighborhoodsTest)); 
  #endif
  #ifdef TEST_2D_ADAPTED_LOCALMINMAX
	add( testCase( &LocalMinTests::test2D_AdaptedLocalminmaxTest)); 
  #endif

	add( testCase( &LocalMinTests::test2D_LocalMinMaxTest)); 
#endif

#ifdef CHECK_NEIGHBORHOODS
	add( testCase( &LocalMinTests::test3D_NeighborhoodCheckTest)); 
#endif


#ifdef TEST_3D_REFERENCE_IMPL
	add( testCase( &LocalMinTests::test3D_ReferenceImplTest)); 
#endif
#ifdef TEST_3D_REFERENCE_IMPL2
	add( testCase( &LocalMinTests::test3D_ReferenceImpl2Test)); 
#endif


#ifdef TEST_CHAINED_ITERS  
	add( testCase( &LocalMinTests::test3D_ChainedItersTest)); 
#endif


#ifdef TEST_LISTBASED
  #ifdef TEST_LISTBASED_CHAINED
	add( testCase( &LocalMinTests::test3D_ListBased_ChainedItersTest));
  #endif

  #ifdef TEST_BGL_COMPAT
    #ifdef TEST_BGL_INDIRECT_PROPMAPS 
      	add( testCase( &LocalMinTests::test3D_BGL_IndirectPropmapsTest));
    #endif
    #ifdef TEST_BGL_DIRECT_PROPMAPS
	// could further split this one:
	add( testCase( &LocalMinTests::test3D_BGL_DirectPropmapsTest));
    #endif
  #endif

#endif
    }
};

int main(int argc, char ** argv) {


      
//     } catch (vigra::StdException &e) {
//       cout << e.what() << endl;
//       return 1;
//     }



#if 1
    GridGraphsTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
#endif


    return 0;
}



/// Local Variables: 
/// c-basic-offset: 4
/// End:
