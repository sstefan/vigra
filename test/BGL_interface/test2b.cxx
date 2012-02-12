// 3D tests:
#define TEST_3D_REFERENCE_IMPL
//#define TEST_3D_REFERENCE_IMPL2 // vigra-side part broken
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
#include <vigra/voxelneighborhood.hxx>
#include <vigra/localminmax.hxx>

#include <vigra/timing.hxx>

#include "empty.hxx"

#ifdef TEST_LISTBASED
// need to include BGL-related bits these before the algorithms:
#include <vigra/multi_gridgraph.hxx>
#include <vigra/multi_gridgraph_coords.hxx>
#endif


// for multi-dim image experiments:
#include <vigra/multi_array.hxx>
#include <vigra/multi_localminmax.hxx>

#ifdef TEST_BGL_COMPAT
#include <boost/graph/graph_concepts.hpp>
#endif


USETICTOC

using namespace std;
using namespace vigra;

const bool verbose = true;
const size_t repetitions = 5;
//const size_t best_of_repetitions = 3;
const size_t best_of_repetitions = 10;


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






struct LocalMinTests : TestImageReader {

#if 1
    // test a variant: 
    //     vertex_descriptor = coord object -> used to access property map.
    //     vertex_iterator = coupled_iterator<0 bd> -> only knows SOI, coord, shape  (an "index"-type object)
    //     adjacency_iterator = [nb-index, nb-ref-list, ref_pt=vertex_iterator]
    void test3D_BGL_DirectPropmapsTest_CoordsVariant()
      { 
	cout << "Testing BGL compatibility on grid graph; Coords variant." << endl;
	TestData_3D td(in, 8);
	using namespace boost;

	typedef image_2D_type::value_type value_type;
	typedef GridGraphView_CoordsDescriptor<3> grid_view_type;
	grid_view_type ggv(td.vol_src.shape());
	//	grid_view_type ggv_dest(td.vol_dest);

	// FIXME: Make both images  properties  of the common underlying grid graph.

	MultiArrayView_property_map<TestData_3D::multiarrayview3d_type> weights(td.vol_src);
	MultiArrayView_property_map<TestData_3D::multiarrayview3d_type> markers(td.vol_dest);

	// access markers/weights "property maps" via corresponding get/set functions
	{
	  cout << "Timing 3D (new, list-based, via BGL (vertex_desc=coord), direct propmaps src/dest): " << endl;
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	      localMinMaxGraph3(ggv,
				weights, // =weights,
				markers, //=markers,
				value_type(255.0), // marker
				value_type(0.0), // threshold
				std::greater<value_type>());
	  TICTOCLOOP_END
	      // save_volume(td.vol_dest,  "slice_wrong_res");
	  cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
	}

#if 1
	// variant: adjacent_vertices(vertex_iterator)
	{
	  cout << "Timing 3D (new, list-based, via BGL (vertex_desc=coord), direct propmaps src/dest, adjacent_vertices(vertex_iterator) variant): " << endl;
	  TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	      localMinMaxGraph4(ggv,
				weights, // =weights,
				markers, //=markers,
				value_type(255.0), // marker
				value_type(0.0), // threshold
				std::greater<value_type>());
	  TICTOCLOOP_END
	  cout << (td.check() ? "Result ok. " : "FAILED! ") << endl;
	}
#endif
      }



#endif 

};




struct GridGraphsTestSuite
: public vigra::test_suite
{
    GridGraphsTestSuite()
    : vigra::test_suite("GridGraphsTestSuite")
    {
	add( testCase( &LocalMinTests::test3D_BGL_DirectPropmapsTest_CoordsVariant));
    }
};

int main(int argc, char ** argv) {


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
