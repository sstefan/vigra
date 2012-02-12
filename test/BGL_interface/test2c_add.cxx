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

// #include <unittest.hxx>
#include <vigra/imageinfo.hxx> // this includes multi_iterator.hxx 
#include <vigra/impex.hxx>
#include <vigra/copyimage.hxx>
#include <vigra/voxelneighborhood.hxx>
#include <vigra/localminmax.hxx>

#include <vigra/timing.hxx>

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


//USETICTOC

using namespace std;
using namespace vigra;

namespace {
const bool verbose = true;
const size_t repetitions = 5;
//const size_t best_of_repetitions = 3;
const size_t best_of_repetitions = 10;

typedef vigra::MultiArray<2, unsigned char> image_2D_type;
typedef MultiArrayShape<2>::type Shp2D;



} // unnamed namespace

struct TestData_3D 
{
  typedef vigra::MultiArray<3, image_2D_type::value_type> multiarray_type;
  typedef multiarray_type::view_type multiarrayview3d_type;
//   typedef vigra::MultiArrayView<3, vigra::BImage::PixelType> multiarrayview3d_type;
  typedef vigra::MultiArrayView<2, image_2D_type::value_type> multiarrayview_type;

  typedef multiarray_type::value_type value_type;

//   const image_2D_type &in_;
//   image_2D_type out;
//   const bool allowAtBorder_;
//   multiarray_type vol_src, vol_dest;
//   multiarray_type::size_type sz;

  bool check(bool verbose=false);
};

struct AddSomeCode {

    typedef image_2D_type::value_type value_type;
    typedef GridGraphView_CoordsDescriptor<3> grid_view_type;

    static void acceleratorCode(TestData_3D &td,
			 grid_view_type &ggv, 
			 boost::MultiArrayView_property_map<TestData_3D::multiarrayview3d_type> &weights,
			 boost::MultiArrayView_property_map<TestData_3D::multiarrayview3d_type> &markers
			 );
};



void AddSomeCode::acceleratorCode(TestData_3D &td,
			 grid_view_type &ggv, 
			 boost::MultiArrayView_property_map<TestData_3D::multiarrayview3d_type> &weights,
			 boost::MultiArrayView_property_map<TestData_3D::multiarrayview3d_type> &markers
			 )
      { 
	cout << "Testing BGL compatibility on grid graph; Coords variant." << endl;
	using namespace boost;

#if 1
	// somehow makes the above test faster... even if only run for 1 iteration
	bool flag = false; // this suffices to compile the below code but not run it, which makes the above run faster
	// const bool flag = false; // this does not suffice, because the block below is optimized away then.
	// variant: adjacent_vertices(vertex_iterator)
	flag = true; 
	if (flag)
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



/// Local Variables: 
/// c-basic-offset: 4
/// End:
