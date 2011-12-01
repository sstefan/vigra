#pragma once
#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/imageinfo.hxx> // this includes multi_iterator.hxx 


// {{{ struct TestData2_3D 

template<class T, int SLICES>
struct TestData2_3D 
{
    typedef T value_type;
    typedef vigra::MultiArray<3, T> multiarray_type;
    typedef vigra::MultiArray<2, T> image_2D_type;
    typedef typename multiarray_type::view_type multiarrayview3d_type;
    typedef vigra::MultiArrayView<2, T> multiarrayview_type;
    
    
    const image_2D_type &in_;
    image_2D_type out;
    multiarray_type vol_src, vol_dest;
    typename multiarray_type::size_type sz;
    
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
	  (vol_src.bind<2>(z)).copy(v_src);
      }
    }
};

// }}}

// {{{ struct TestData3_nD

template<class T, int N>
struct TestData3_nD
{
    typedef T value_type;
    typedef vigra::MultiArray<N, T> multiarray_type;
    typedef typename multiarray_type::view_type view_type;
    
    typename multiarray_type::size_type sz;

    multiarray_type vol_src;
    //    multiarray_type  vol_dest;

	view_type vol_view;

    typedef vigra::TinyVector<int,N> IntVec;
    typedef vigra::TinyVector<double,N> DoubleVec;

    enum { dim = N };
    //    enum { SX = 100, SY=100, SZ = 20 }; // dimensions

    TestData3_nD() 
	: sz(20), vol_src(sz, 0.0), vol_view(vol_src)
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
void save_volume(const vigra::MultiArrayView<3, T> &view, std::string fnbase = std::string("slice")) 
{
    for (size_t z=0; z<view.shape()[2]; ++z) {
	std::stringstream fn;
	fn << fnbase << z << ".png";
	//	exportImage(srcImageRange(view.bind<2>(z)), vigra::ImageExportInfo(fn.str().c_str())); // clang doesn't understand this
 	exportImage(srcImageRange(view.bindOuter(z)), vigra::ImageExportInfo(fn.str().c_str()));
    }
}
