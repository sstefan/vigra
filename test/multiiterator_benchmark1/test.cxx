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
// #define NONGRID // non-grid-structured graph test case

#define TEST_WATERSHEDS_3D_REF

#include <typeinfo>

#include <vigra/stdimage.hxx>
#include <vigra/basicimageview.hxx>

typedef vigra::BasicImage<double> ImageType;
typedef vigra::BasicImageView<double> ImageViewType;

#include "../test/include/unittest.hxx"

#include <iostream>
#include <boost/program_options.hpp>

#ifdef NONGRID
#include "graphtest.hxx"
#endif


#include <vigra/imageinfo.hxx> // this includes multi_iterator.hxx 
#include <vigra/impex.hxx>
#include <vigra/copyimage.hxx>
// #include <vigra/edgedetection.hxx>
//#include <vigra/localminmax.hxx>
#include <vigra/voxelneighborhood.hxx>
#include <vigra/localminmax.hxx>

#include <vigra/timing.hxx>

// for multi-dim image experiments:
#include <vigra/multi_array.hxx>
#include <vigra/multi_localminmax.hxx>

#ifdef TEST_LISTBASED
#include <vigra/multi_gridgraph.hxx>
#endif


#ifdef TEST_BGL_COMPAT
#include <boost/graph/graph_concepts.hpp>
#endif


// Watersheds:
#ifdef TEST_WATERSHEDS_3D_REF
#include <vigra/multi_pointoperators.hxx>
//#include <vigra/watersheds3d.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/seededregiongrowing3d.hxx>
#endif

USETICTOC

using namespace std;
using namespace vigra;


typedef vigra::MultiArray<2, unsigned char> image_2D_type;
typedef MultiArrayShape<2>::type Shp2D;


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

// {{{ struct TestData2_3D 

template<int SLICES>
struct TestData2_3D 
{
    typedef vigra::MultiArray<3, image_2D_type::value_type> multiarray_type;
    typedef multiarray_type::view_type multiarrayview3d_type;
    typedef vigra::MultiArrayView<2, image_2D_type::value_type> multiarrayview_type;
    
    typedef multiarray_type::value_type value_type;
    
    const image_2D_type &in_;
    image_2D_type out;
    multiarray_type vol_src, vol_dest;
    multiarray_type::size_type sz;
    
    TestData2_3D(const image_2D_type &in)
	: in_(in), out(in.shape())
    {
      multiarrayview_type v_src(in);
      
      sz[2]=SLICES;
      sz[0]=in.shape()[0];
      sz[1]=in.shape()[1];

      vol_src = multiarray_type(sz, 0.0);
      vol_dest = multiarray_type(sz, 0.0);
      for (size_t z=0; z<SLICES; z+=1) {
	  vol_src.bind<2>(z).copy(v_src);
      }
    }
};

// }}}

// {{{ struct TestData3_nD

template<int N>
struct TestData3_nD
{
    typedef vigra::MultiArray<N, image_2D_type::value_type> multiarray_type;
    typedef typename multiarray_type::view_type view_type;
    typedef typename multiarray_type::value_type value_type;
    
    typename multiarray_type::size_type sz;

    multiarray_type vol_src;
    //    multiarray_type  vol_dest;

    typedef vigra::TinyVector<int,N> IntVec;
    typedef vigra::TinyVector<double,N> DoubleVec;

    enum { dim = N };
    //    enum { SX = 100, SY=100, SZ = 20 }; // dimensions

    TestData3_nD() 
	: sz(20), vol_src(sz, 0.0)
    {
	// use a scan-order iterator to fill volume
	// with the minimum of the Euclidean distance to a set of test seed points:
	// (located at 1/3, 2/3 along the hyperdiagonal(s))

	DoubleVec pt1 = vol_src.shape()*(1.0/3.0);
	DoubleVec pt2 = vol_src.shape()*(2.0/3.0);

	typename view_type::iterator it=view_type(vol_src).begin();
	typename view_type::iterator itend=view_type(vol_src).end();

	for (;it != itend; ++it) {
	    IntVec temp = it.point();
	    double y1 = (temp-pt1).squaredMagnitude();
	    double y2 = (temp-pt2).squaredMagnitude();
	    *it = y1 < y2 ? y1 : y2;
	}
    }
};

// }}}


template<class T>
void save_volume(const vigra::MultiArrayView<3, T> &view, string fnbase = string("slice")) 
{
    for (size_t z=0; z<view.shape()[2]; ++z) {
	stringstream fn;
	fn << fnbase << z << ".png";
	//	exportImage(srcImageRange(view.bind<2>(z)), vigra::ImageExportInfo(fn.str().c_str())); // clang doesn't understand this
	exportImage(srcImageRange(view.bindOuter(z)), vigra::ImageExportInfo(fn.str().c_str()));
    }
}




#if 1 
struct GridGraphsTestSuite
: public vigra::test_suite
{
    GridGraphsTestSuite()
    : vigra::test_suite("GridGraphsTestSuite")
    {
	//        add( testCase( &Watersheds3dTest::testDistanceVolumesSix));
    }
};
#endif

int main(int argc, char ** argv) {
    const bool verbose = true;
    const size_t repetitions = 1;
    const size_t best_of_repetitions = 3;

    string input("test.png");
    string output("out.png");

    try {
      // read image using vigra
      if (!isImage(input.c_str())) {
	cerr << "Input image not readable." << endl;
	return 1;
      }
      cout << "Reading image " << input << endl;
      vigra::ImageImportInfo info(input.c_str());

      Shp2D imsize(info.size());
      //      Shp2D imsize(info.width(), info.height());
      cout << "INPUT IMAGE SIZE: " << imsize << endl;

      image_2D_type in(imsize);
      image_2D_type out(imsize);

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

      // init output
      out *= 0;

      cout << "Running "<<repetitions <<  "x" << best_of_repetitions << " repetitions:" << endl;
      TIC;
      for (int i=0; i<repetitions; i++)
      {
	// local min/max extraction
	vigra::localMaxima(
			   srcImageRange(in), 
			   destImage(out),
			   //			   vigra::LocalMinmaxOptions().neighborhood(4).allowAtBorder().threshold(0).markWith(255)
			   //vigra::LocalMinmaxOptions().neighborhood(8).allowAtBorder().threshold(0).markWith(255)
			   //			   vigra::LocalMinmaxOptions().neighborhood(8).threshold(0).markWith(255)
			   vigra::LocalMinmaxOptions().neighborhood(4).threshold(0).markWith(255)
			   );
      }
      double time=TOCN;
      if (verbose) 
	cout << "Timing: " << repetitions << " repetitions took " << time << " ms." << endl;

      // again using new timing code:
      TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	vigra::localMaxima(
			   srcImageRange(in), 
			   destImage(out),
			   //			   vigra::LocalMinmaxOptions().neighborhood(4).allowAtBorder().threshold(0).markWith(255)
			   //vigra::LocalMinmaxOptions().neighborhood(8).allowAtBorder().threshold(0).markWith(255)
			   //			   vigra::LocalMinmaxOptions().neighborhood(8).threshold(0).markWith(255)
			   vigra::LocalMinmaxOptions().neighborhood(4).threshold(0).markWith(255)
			   );
      TICTOCLOOP_END

      cout << "Writing result " << output << endl;
      exportImage(srcImageRange(out), vigra::ImageExportInfo("out0.png"));

      out *= 0;

      // new code (testing):
      // get a MultiArrayView on the input image
      typedef vigra::MultiArrayView<2, image_2D_type::value_type> multiarrayview_type;
      multiarrayview_type v_src = multiarrayview_type(in);
      multiarrayview_type v_dest = multiarrayview_type(out);
      cout << "SOURCE VIEW SIZE:" << v_src.shape() << endl;

#ifdef TEST_2D      
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
#endif


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


      cout << "Writing result " << output << endl;
      exportImage(srcImageRange(out), vigra::ImageExportInfo(output.c_str()));

#endif // TEST_2D


#ifdef CHECK_NEIGHBORHOODS
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
      // this code is commented out in Vigra for an unknown reason.
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


#ifdef TEST_WATERSHEDS_3D_REF
      // Watershed Tests:
      // (adapted from Vigra example for 2D)
      {
	cout << "Testing Watershed code:" << endl;
	  // FIXME: Use a different (synthetic!) test image here, 
	  // which actually has local 3D minima!
	//TestData2_3D<5> td(in);
	TestData3_nD<3> td;
	save_volume(td.vol_src, "testvol");
	  // a labeling image of same size
	  typedef vigra::MultiArray<3, unsigned int> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  // TestData_3D::multiarrayview3d_type gradMag(td.vol_src);
	  TestData3_nD<3>::view_type gradMag(td.vol_src);
	  unsigned int max_region_label;

	  cout << "Timing 3D SRG(reference impl.): " << endl;
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
		  // compute seeds beforehand (use connected components of all pixels 
		  // where the gradient  is below 4.0)
		  max_region_label = 
	      generateWatershedSeeds3D(srcMultiArrayRange(gradMag), 
				       destMultiArray(labeling),
				       NeighborCode3DTwentySix(),
				       SeedOptions().levelSets(4.0));
// 		      generateWatershedSeeds3D(gradMag, labelarray_type::view_type(labeling),
// 					       NeighborCode3DTwentySix(),
// 					       SeedOptions().levelSets(4.0));
		  cout << "Max region label: " << max_region_label << endl;
		  save_volume(labeling, "labeling");

		  // quantize the gradient image to 256 gray levels
		  typedef vigra::MultiArray<3, unsigned char> quantized_type;
		  quantized_type gradMag256(td.vol_src.shape(), static_cast<unsigned char>(0));

		  vigra::FindMinMax<float> minmax; 
		  inspectMultiArray(srcMultiArrayRange(gradMag), minmax); // find original range
		  transformMultiArray(srcMultiArrayRange(gradMag), destMultiArray(gradMag256),
				      linearRangeMapping(minmax, 0, 255));
		  save_volume(gradMag256, "gradMag256");

		  // create a statistics functor for region growing
		  detail::WatershedStatistics<quantized_type::value_type, unsigned int> regionstats;

		  seededRegionGrowing3D(srcMultiArrayRange(gradMag256), srcMultiArray(labeling),
					destMultiArray(labeling), regionstats, KeepContours);
	      TICTOCLOOP_END
	      cout << "  (Including file I/O)" << endl;
	  // cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
        
	  // FIXME: The 2D variant of watershedsRegionGrowing includes the seed computation???
	  // (only when option "mini" is not set to "Unspecified".)

	  // 	  // call the turbo algorithm with 256 bins, using 8-neighborhood
	  // 	  watershedsRegionGrowing(srcImageRange(gradMag256), destImage(labeling),
	  // 				  WatershedOptions().turboAlgorithm(256));
	  cout << "Max region label after region growing: " << max_region_label << endl;
	  save_volume(labeling, "labeling_res");


	  // Again, using the watershedsRegionGrowing_nD code
	  // (not really nD yet, but soon...)
#if 1 
	  cout << "watershedsRegionGrowing_nD():" << endl;
	  max_region_label = 
	      watershedsRegionGrowing_nD(srcMultiArrayRange(gradMag),
					 destMultiArray(labeling),
 					 NeighborCode3DTwentySix(), // FIXME generalize!
					 WatershedOptions()); // FIXME: Test different options
	  cout << "Max region label after watershedsRegionGrowing_nD: " << max_region_label << endl;
#endif
      }
#endif


#define TEST_WATERSHEDS_3D
#ifdef TEST_WATERSHEDS_3D
      {
	  const int dim=3;
	cout << "Testing Watershed code:" << endl;
	  // FIXME: Use a different (synthetic!) test image here, 
	  // which actually has local 3D minima!
	TestData3_nD<dim> td;
	//	save_volume(td.vol_src, "testvol");
	  // a labeling image of same size
	  typedef vigra::MultiArray<dim, unsigned int> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  TestData3_nD<dim>::view_type gradMag(td.vol_src);
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


#define TEST_WATERSHEDS_4D
#ifdef TEST_WATERSHEDS_4D
      {
	  const int dim=4;
	cout << "Testing Watershed code:" << endl;
	  // FIXME: Use a different (synthetic!) test image here, 
	  // which actually has local 3D minima!
	TestData3_nD<dim> td;
	//	save_volume(td.vol_src, "testvol");
	  // a labeling image of same size
	  typedef vigra::MultiArray<dim, unsigned int> labelarray_type;
	  labelarray_type labeling(td.sz, 0.0);
        
	  // compute a gradient magnitude volume for input:
	  // TODO; for now, copy image itself
	  TestData3_nD<dim>::view_type gradMag(td.vol_src);
	  unsigned int max_region_label;

	  cout << "watershedsRegionGrowing_nD():" << endl;
// 	  max_region_label = 
// 	      watershedsRegionGrowing_nD(srcMultiArrayRange(gradMag),
// 					 destMultiArray(labeling),
//  					 NeighborCode3DTwentySix(), // FIXME generalize!
// 					 WatershedOptions()); // FIXME: Test different options
	  max_region_label = 
	      watershedsRegionGrowing_nD(gradMag,
					 labeling,
 					 NeighborCode3DTwentySix(), // FIXME generalize!
					 WatershedOptions()); // FIXME: Test different options
	  cout << "Max region label after watershedsRegionGrowing_nD: " << max_region_label << endl;
      }
#endif

      
    } catch (vigra::StdException &e) {
      cout << e.what() << endl;
      return 1;
    }

    // {{{ NONGRID test

#ifdef NONGRID
    if (verbose) 
      cout << "Testing Non-Grid-Graph:" << endl;

    { 
      // test a non-grid graph (via BGL)
      using namespace boost;

      typedef StS::TestGraph::Graph Graph;
      StS::TestGraph TG;
      TG.test();

      if (verbose) 
	cout << "Testing Non-Grid-Graph LocalMax:" << endl;

      // create an associated property map taking bools to indicate 
      //   if node is locally maximal:

      // to apply the same localminmax algorithm,
      // require G to be a model of AdjacencyGraph.
      // Furthermore, we need a method to associate the boolean label with each node,
      // i.e. a property map over the nodes of the graph.

      // external properties:
      vector<bool> markers_vec(num_vertices(TG.G));
      vector<double> &weights_vec(TG.weights_); 
      // TODO: Also test internal property variant for weights!

      // define corresponding property maps:
      typedef property_map<Graph, vertex_index_t>::type IndexMap;
      IndexMap index = get(vertex_index, TG.G);

      // iterator_property_map is responsible for mapping nodes to keys
      iterator_property_map<vector<bool>::iterator, IndexMap>
	markers(markers_vec.begin(), index);
      iterator_property_map<vector<double>::iterator, IndexMap>
	weights(weights_vec.begin(), index);

      localMinMaxGraph(TG.G,
		       weights,
		       markers,
		       true, // marker
		       (double)0.0, // threshold
		       std::greater<double>());

      // Output result
      typedef boost::graph_traits<StS::TestGraph::Graph>::vertex_iterator 
	graph_scanner;
      graph_scanner it, itend;
      cout << "Graph LocalMax Result:" << endl;
      for (tie(it,itend)=vertices(TG.G); it!=itend; ++it) {
	cout << TG.names_[*it];
      }
      cout << endl;
      for (tie(it,itend)=vertices(TG.G); it!=itend; ++it) {
	cout << (markers[*it] ? "M" : ".");
      }
      cout << endl;
      bool expected_markers[] = { true, false, false, false, true, false, true, 
				  false, false, false, false, false, true };
      bool check_failed = false;
      int nr=0;
      for (tie(it,itend)=vertices(TG.G); it!=itend; ++it, ++nr) {
	check_failed |= (expected_markers[nr] != markers[*it]);
      }
      cout << (check_failed ? "FAILED." : "OK") << endl;
    }
#endif //  NONGRID

    // }}}




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
