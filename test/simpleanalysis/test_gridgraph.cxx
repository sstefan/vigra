/************************************************************************/
/*                                                                      */
/*    Copyright 2004-2012 by Stefan Schmidt, Ullrich Koethe             */
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

#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/distancetransform.hxx"
#include "vigra/cornerdetection.hxx"
#include "vigra/symmetry.hxx"
#include "vigra/noise_normalization.hxx"
#include "vigra/affinegeometry.hxx"
#include "vigra/affine_registration.hxx"
#include "vigra/impex.hxx"

#include <vigra/multi_array.hxx>
#include <vigra/multi_gridgraph_coords.hxx>
#include <vigra/multi_localminmax.hxx>
#include <vigra/multi_labelgraph.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/multi_seededregiongrowing.hxx>

#ifdef HasFFTW3
# include "vigra/slanted_edge_mtf.hxx"
#endif

using namespace vigra;
using namespace vigragraph;

struct LabelingTest
{
    typedef MultiArray<2, double> Image;
    typedef Image::difference_type Shape;
    typedef MultiArrayView<2, Image::value_type> ViewType;

    typedef GridGraphView_CoordsDescriptor<2> GridGraph;
    typedef MultiArrayView_property_map<ViewType> VertexPropMap;

    LabelingTest()
	: img1(Shape(5,5)), img2(Shape(5,5)), img3(Shape(9,5)), img4(Shape(11,11))
    {
        static const double in1[] = { 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0};

	ViewType::iterator i = static_cast<ViewType>(img1).begin();
        ViewType::iterator end = static_cast<ViewType>(img1).end();
        const double * p = in1;

        for(; i != end; ++i, ++p)
        {
            *i = *p;
        }

        static const double in2[] = { 0.0, 1.0, 0.0, 1.0, 0.0,
                                      1.0, 0.0, 1.0, 0.0, 1.0,
                                      0.0, 1.0, 0.0, 1.0, 0.0,
                                      1.0, 0.0, 1.0, 0.0, 1.0,
                                      0.0, 1.0, 0.0, 1.0, 0.0};

        i = static_cast<ViewType>(img2).begin();
        end = static_cast<ViewType>(img2).end();
        p = in2;

        for(; i != end; ++i, ++p)
        {
            *i = *p;
        }

        static const int in3[] = {
            0, 1, 0, 1, 0, 1, 0, 1, 0,
            0, 1, 0, 1, 0, 1, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 0, 1, 0,
            0, 1, 0, 0, 0, 1, 1, 1, 0,
            0, 0, 0, 1, 1, 1, 0, 0, 0
        };

        i = static_cast<ViewType>(img3).begin();
        end = static_cast<ViewType>(img3).end();
        const int * p1 = in3;

        for(; i != end; ++i, ++p1)
        {
            *i = *p1;
        }

        static const int spiral[] = {
            1,1,1,1,1,1,1,1,1,1,1,
            1,2,2,2,2,2,2,2,2,2,1,
            1,2,1,1,1,1,1,1,1,2,1,
            1,2,1,2,2,2,2,2,1,2,1,
            1,2,1,2,1,1,1,2,1,2,1,
            1,2,1,2,1,2,1,2,1,2,1,
            1,2,1,2,1,2,2,2,1,2,1,
            1,2,1,2,1,1,1,1,1,2,1,
            1,2,1,2,2,2,2,2,2,2,1,
            1,2,1,1,1,1,1,1,1,1,1,
            1,2,2,2,2,2,2,2,2,2,2
        };

        i = static_cast<ViewType>(img4).begin();
        end = static_cast<ViewType>(img4).end();
        p1 = spiral;

        for(; i != end; ++i, ++p1)
        {
            *i = *p1;
        }
    }

    void labelingFourTest1()
    {
        Image res(img1);
	GridGraph graph(img1.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	int count = labelGraph(graph, VertexPropMap(static_cast<ViewType>(img1)), out);

        should(2 == count);

        ViewType::iterator i1 = static_cast<ViewType>(img1).begin();
        ViewType::iterator i1end = static_cast<ViewType>(img1).end();
        ViewType::iterator i2 = static_cast<ViewType>(res).begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == *i2 - 1.0);
        }
    }

    void labelingFourTest2()
    {
        Image res(img3);
	GridGraph graph(img3.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	int count = labelGraph(graph, VertexPropMap(static_cast<ViewType>(img3)), out);

        should(6 == count);

        static const int target[] = {
            1, 2, 1, 3, 1, 4, 1, 5, 1,
            1, 2, 1, 3, 1, 4, 1, 1, 1,
            1, 2, 1, 3, 1, 1, 1, 6, 1,
            1, 2, 1, 1, 1, 6, 6, 6, 1,
            1, 1, 1, 6, 6, 6, 1, 1, 1
        };

        ViewType::iterator i = static_cast<ViewType>(res).begin();
        ViewType::iterator iend = static_cast<ViewType>(res).end();
        const int * p = target;

        for(; i != iend; ++i, ++p)
        {
            should(*i == *p);
        }
    }

    void labelingFourTest3()
    {
        Image res(img4.shape());
	GridGraph graph(img4.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	should(2 == labelGraph(graph, VertexPropMap(static_cast<ViewType>(img4)), out));

        ViewType::iterator i = static_cast<ViewType>(res).begin();
        ViewType::iterator iend = static_cast<ViewType>(res).end();
        ViewType::iterator id = static_cast<ViewType>(img4).begin();

        for(; i != iend; ++i, ++id)
        {
            should(*i == *id);
        }
    }

    void labelingFourTest4()
    {
        static const int data[] = {
            1,1,1,1,1,1,1,1,1,2,
            2,1,1,1,1,1,1,2,1,2,
            2,1,1,1,1,2,1,2,1,2,
            2,1,1,2,1,2,1,2,1,2,
            2,1,2,2,2,2,2,2,2,2,
            2,2,2,3,3,3,3,3,3,3

        };

        int w=10;
        int h=6;
        Image img(Shape(w,h)), res(Shape(w,h));

        std::copy(data, data+w*h, static_cast<ViewType>(img).begin());

	GridGraph graph(img.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	should(3 == labelGraph(graph, VertexPropMap(static_cast<ViewType>(img)), out));

        ViewType::iterator i = static_cast<ViewType>(res).begin();
        ViewType::iterator iend = static_cast<ViewType>(res).end();
	ViewType::iterator id = static_cast<ViewType>(img).begin();

        for(int c=0; i != iend; ++i, ++id, ++c)
        {
            should(*i == *id);
        }
    }

    void labelingToCrackEdgeTest()
    {
        Image tmp(img1);
        Image res(Shape(9, 9));

	GridGraph graph(tmp.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(tmp));
	should(2 == labelGraph(graph, VertexPropMap(static_cast<ViewType>(img1)), out));

        regionImageToCrackEdgeImage(srcImageRange(tmp), destImage(res), 0.0);

        static const double desired[] = {
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
               1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,
               1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

        const double * i1 = desired;
        const double * i1end = i1 + 81;
        ViewType::iterator i2 = static_cast<ViewType>(res).begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == *i2);
        }
    }

    void labelingEightTest1()
    {
        Image res(img2);

	GridGraph graph(img2.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	should(2 == labelGraph(graph, VertexPropMap(static_cast<ViewType>(img2)), out));

        ViewType::iterator i1 = static_cast<ViewType>(img2).begin();
        ViewType::iterator i1end = static_cast<ViewType>(img2).end();
        ViewType::iterator i2 = static_cast<ViewType>(res).begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == *i2 - 1.0);
        }
    }

    void labelingEightTest2()
    {
        static const int data[] = {
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,2,3,3,2,1,1,1,
            1,1,2,3,3,2,4,2,1,1,
            1,2,3,3,2,4,4,4,2,1,
            1,2,3,2,4,4,2,3,2,1,
            1,2,3,3,2,2,3,3,2,1,
            1,1,2,3,3,3,3,2,1,1,
            1,1,1,2,2,2,2,1,1,1,
            1,1,1,1,1,1,1,1,1,1
        };

        Image img(Shape(10,9)), res(Shape(10,9));

        std::copy(data, data+90, static_cast<ViewType>(img).begin());

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	should(4 == labelGraph(graph, VertexPropMap(static_cast<ViewType>(img)), out));

        ViewType::iterator i = static_cast<ViewType>(res).begin();
        ViewType::iterator iend = static_cast<ViewType>(res).end();
        const int * p = data;

        for(; i != iend; ++i, ++p)
        {
            should(*i == *p);
        }
    }

    void labelingFourWithBackgroundTest1()
    {
        Image res(img1);
	res.init(0);

	GridGraph graph(img1.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	should(1 == labelGraphWithBackground(graph, 
					     VertexPropMap(static_cast<ViewType>(img1)), 
					     out,
					     0.0));

        ViewType::iterator i1 = static_cast<ViewType>(img1).begin();
        ViewType::iterator i1end = static_cast<ViewType>(img1).end();
        ViewType::iterator i2 = static_cast<ViewType>(res).begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == *i2);
        }
    }

    void labelingFourWithBackgroundTest2()
    {
        Image res(img4);

	GridGraph graph(img4.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));

	should(1 == labelGraphWithBackground(graph, 
					     VertexPropMap(static_cast<ViewType>(img4)), 
					     out,
					     2));

        ViewType::iterator i = static_cast<ViewType>(res).begin();
        ViewType::iterator iend = static_cast<ViewType>(res).end();
        ViewType::iterator id = static_cast<ViewType>(img4).begin();

        for(; i != iend; ++i, ++id)
        {
            should(*i == (2-*id));
        }
    }

    void labelingEightWithBackgroundTest()
    {
        Image res(img2);

	GridGraph graph(img2.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	should(1 == labelGraphWithBackground(graph, 
					     VertexPropMap(static_cast<ViewType>(img2)), 
					     out,
					     0.0));

        ViewType::iterator i1 = static_cast<ViewType>(img2).begin();
        ViewType::iterator i1end = static_cast<ViewType>(img2).end();
        ViewType::iterator i2 = static_cast<ViewType>(res).begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == *i2);
        }
    }

    Image img1, img2, img3, img4;
};



template <class T>
struct EqualWithToleranceFunctor
{
    EqualWithToleranceFunctor(T tolerance = 2.0 * NumericTraits<T>::epsilon())
    : t(tolerance)
    {}
    
    bool operator()(T l, T r) const
    {
        return abs(l-r) <= t;
    }
    
    T t;
};

struct LocalMinMaxTest
{
    typedef MultiArray<2, double> Image;
    typedef Image::difference_type Shape;
    typedef MultiArrayView<2, Image::value_type> ViewType;
    typedef MultiArray<2, Int32> IImage;
    typedef MultiArrayView<2, IImage::value_type> ViewTypeI;

    typedef GridGraphView_CoordsDescriptor<2> GridGraph;
    typedef MultiArrayView_property_map<ViewType> VertexPropMap;
    typedef MultiArrayView_property_map<ViewTypeI> VertexPropMapI;

    typedef vigra::MultiArray<3,double> Volume;
    typedef vigra::MultiArray<3,Int32> IVolume;
    typedef MultiArrayShape<3>::type Shp3D;
    typedef MultiArrayView<3, Volume::value_type> ViewType3;
    typedef MultiArrayView<3, IVolume::value_type> ViewType3I;

    typedef GridGraphView_CoordsDescriptor<3> GridGraph3;
    typedef MultiArrayView_property_map<ViewType3> VertexPropMap3;
    typedef MultiArrayView_property_map<ViewType3I> VertexPropMap3I;

    LocalMinMaxTest()
	: img(Shape(9,9)), vol()
    {
        static const double in[] = {
            0.2,  0.1,  0.1,  0.3,  0.5,  0.3,  0.0,  0.0, -0.1,
            0.0, -0.1,  0.1,  0.0,  1.0,  0.0,  0.3,  0.0,  0.0,
            0.0,  0.5,  2.0,  0.0,  2.0,  2.0,  2.0,  0.0, -1.1,
            0.1,  0.0,  1.0,  1.5,  1.0,  1.0,  0.0,  0.0,  0.0,
            0.0,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
            0.0,  0.0,  0.0,  0.0, -1.0, -1.5, -1.0,  0.0, -0.3,
            0.0,  0.0, -2.0, -2.0, -2.0,  0.0, -2.0, -0.5,  0.0,
            1.0,  0.0,  0.0,  0.0, -1.0,  0.0, -0.1,  0.1,  0.0,
            0.0,  0.0,  0.0,  0.0, -0.5, -0.3, -0.1, -0.1,  0.0};

        ViewType::iterator i = static_cast<ViewType>(img).begin();
        ViewType::iterator end = static_cast<ViewType>(img).end();
        const double * p = in;

        for(; i != end; ++i, ++p)
        {
            *i = *p;
        }

        //prepare the multiarray
        vol.reshape(Shp3D(10,20,50),0);

        vol(1,1,1)=10;
        vol(5,5,5)=350;

        vol(8,3,5)=9; //plateau
        vol(8,4,5)=9;
        vol(8,5,5)=9;

        vol(6,7,7)=-0.5;
        vol(7,7,7)=-1;
        vol(7,1,15)=-100;
        vol(7,1,19)=-20;

        vol(3,15,26)=-1; //plateau
        vol(3,15,27)=-1;
        vol(3,15,28)=-1;
        vol(3,16,26)=-1;

	vol(9,18,35)=-100; //on the border is NOT skipped
	vol(0,1,49)=100; //on the border is NOT skipped

    }


    void localMinimum3DTest()
    {
        Volume res(vol);
        res.init(0);

	GridGraph3 graph(vol.shape(), DirectNeighborhood);
	VertexPropMap3 out(static_cast<ViewType3>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap3(static_cast<ViewType3>(vol)), 
			       out,
			       Volume::value_type(1), // marker
			       NumericTraits<Volume::value_type>::max(), // threshold
			       std::less<Volume::value_type>()
			       );
	// Difference to this implementation: The above also detects border extrema
        // localMinima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(7,7,7)=1;
        desired(7,1,15)=1;
        desired(7,1,19)=1;
	desired(9,18,35)=1; // also detect on border

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x) {
                    shouldEqual(res(x,y,z), desired(x,y,z));
		}

    }

    void extendedLocalMinimum3DTest()
    {
        Volume res(vol);
        res.init(0);

        IVolume tmp(vol);
	tmp.init(0);

	GridGraph3 graph(vol.shape(), DirectNeighborhood);
	VertexPropMap3 out(static_cast<ViewType3>(res));
	VertexPropMap3I tmpview(static_cast<ViewType3I>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap3(static_cast<ViewType3>(vol)), 
				 out,
				 tmpview,
				 Volume::value_type(1), // marker
				 NumericTraits<Volume::value_type>::max(), // threshold
				 std::less<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

	//        extendedLocalMinima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(7,7,7)=1;
        desired(7,1,15)=1;
        desired(7,1,19)=1;

        desired(3,15,26)=1; //plateau
        desired(3,15,27)=1;
        desired(3,15,28)=1;
        desired(3,16,26)=1;

	desired(9,18,35)=1; // also detect on border

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

    }
    
        void extendedLocalMinimum3DTest2()
    {
        Volume res(vol);
        res.init(0);

        IVolume tmp(vol);
	tmp.init(0);

	GridGraph3 graph(vol.shape(), IndirectNeighborhood);
	VertexPropMap3 out(static_cast<ViewType3>(res));
	VertexPropMap3I tmpview(static_cast<ViewType3I>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap3(static_cast<ViewType3>(vol)), 
				 out,
				 tmpview,
				 Volume::value_type(1), // marker
				 NumericTraits<Volume::value_type>::max(), // threshold
				 std::less<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

	//        extendedLocalMinima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DTwentySix());

        Volume desired(vol);
        desired.init(0);

        desired(7,7,7)=1;
        desired(7,1,15)=1;
        desired(7,1,19)=1;

        desired(3,15,26)=1; //plateau
        desired(3,15,27)=1;
        desired(3,15,28)=1;
        desired(3,16,26)=1;

	desired(9,18,35)=1; // also detect on border

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

    }

    void localMaximum3DTest()
    {
        Volume res(vol);
        res.init(0);

	GridGraph3 graph(vol.shape(), DirectNeighborhood);
	VertexPropMap3 out(static_cast<ViewType3>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap3(static_cast<ViewType3>(vol)), 
			       out,
			       Volume::value_type(1), // marker
			       NumericTraits<Volume::value_type>::min(), // threshold
			       std::greater<Volume::value_type>()
			       );

	//        localMaxima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(1,1,1)=1;
        desired(5,5,5)=1;
	desired(0,1,49)=1; // also detect on border

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));
    }

    void extendedLocalMaximum3DTest()
    {
        Volume res(vol);
        res.init(0);

        IVolume tmp(vol);
	tmp.init(0);

	GridGraph3 graph(vol.shape(), DirectNeighborhood);
	VertexPropMap3 out(static_cast<ViewType3>(res));
	VertexPropMap3I tmpview(static_cast<ViewType3I>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap3(static_cast<ViewType3>(vol)), 
				 out,
				 tmpview,
				 Volume::value_type(1), // marker
				 NumericTraits<Volume::value_type>::min(), // threshold
				 std::greater<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

	//      extendedLocalMaxima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(1,1,1)=1;
        desired(5,5,5)=1;

        desired(8,3,5)=1;
        desired(8,4,5)=1;
        desired(8,5,5)=1;

	desired(0,1,49)=1; // also detect on border

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));
    }
    
        void extendedLocalMaximum3DTest2()
    {
        Volume res(vol);
        res.init(0);
        IVolume tmp(vol);
	tmp.init(0);

	GridGraph3 graph(vol.shape(), IndirectNeighborhood);
	VertexPropMap3 out(static_cast<ViewType3>(res));
	VertexPropMap3I tmpview(static_cast<ViewType3I>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap3(static_cast<ViewType3>(vol)), 
				 out,
				 tmpview,
				 Volume::value_type(1), // marker
				 NumericTraits<Volume::value_type>::min(), // threshold
				 std::greater<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

        // extendedLocalMaxima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DTwentySix());

        Volume desired(vol);
        desired.init(0);

        desired(1,1,1)=1;
        desired(5,5,5)=1;

        desired(8,3,5)=1;
        desired(8,4,5)=1;
        desired(8,5,5)=1;

	desired(0,1,49)=1; // also detect on border

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));
    }

    void localMinimumTest()
    {
        Image res(img);
        res.init(0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap(static_cast<ViewType>(img)), 
			       out,
			       Image::value_type(1), // marker
			       NumericTraits<Image::value_type>::max(), // threshold
			       std::less<Volume::value_type>()
			       );

	//        localMinima(srcImageRange(img), destImage(res));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMinima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowAtBorder());
	// only allowAtBorder is implemented so far.

        desired[8] = 1.0;
        desired[26] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMinimum4Test()
    {
        Image res(img);
        res.init(0);

	GridGraph graph(img.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap(static_cast<ViewType>(img)), 
			       out,
			       Image::value_type(1), // marker
			       NumericTraits<Image::value_type>::max(), // threshold
			       std::less<Volume::value_type>()
			       );

	//        localMinima(srcImageRange(img), destImage(res), LocalMinmaxOptions().neighborhood(4));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMinima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().neighborhood(4).allowAtBorder());
        desired[8] = 1.0;
        desired[26] = 1.0;
        desired[53] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMinimumTestThr()
    {
        Image res(img);
        res.init(0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap(static_cast<ViewType>(img)), 
			       out,
			       Image::value_type(1.0), // marker
			       Image::value_type(-1.0), // threshold
			       std::less<Volume::value_type>()
			       );

//         localMinima(srcImageRange(img), destImage(res),
//                     LocalMinmaxOptions().neighborhood(8).markWith(1.0).threshold(-1.0));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMinima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().neighborhood(8).threshold(-1.0).allowAtBorder());
        desired[26] = 1.0;
	shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMaximumTest()
    {
        Image res(img);
        res.init(0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap(static_cast<ViewType>(img)), 
			       out,
			       Image::value_type(1), // marker
			       NumericTraits<Image::value_type>::min(), // threshold
			       std::greater<Volume::value_type>()
			       );

//         localMaxima(srcImageRange(img), destImage(res));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowAtBorder());
        desired[0] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMaximum4Test()
    {
        Image res(img);
        res.init(0);

	GridGraph graph(img.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap(static_cast<ViewType>(img)), 
			       out,
			       Image::value_type(1), // marker
			       NumericTraits<Image::value_type>::min(), // threshold
			       std::greater<Volume::value_type>()
			       );

//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().neighborhood(4));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().neighborhood(4).allowAtBorder());
        desired[0] = 1.0;
        desired[27] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMaximumTestThr()
    {
        Image res(img);
        res.init(0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	// FIXME: provide more convenient interface!
	localMinMaxGraph3boost(graph, 
			       VertexPropMap(static_cast<ViewType>(img)), 
			       out,
			       Image::value_type(1.0), // marker
			       Image::value_type(0.2), // threshold
			       std::greater<Volume::value_type>()
			       );
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().markWith(1.0).neighborhood(8).threshold(0.2));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowAtBorder().threshold(0.2));
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void extendedLocalMinimumTest()
    {
        Image res(img);
        res.init(0);

        IImage tmp(img);
	tmp.init(0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap(static_cast<ViewType>(img)), 
				 out,
				 tmpview,
				 Image::value_type(1.0), // marker
				 NumericTraits<Volume::value_type>::max(), // threshold
				 std::less<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

//         extendedLocalMinima(srcImageRange(img), destImage(res), 1.0);

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMinima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowPlateaus());
//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMinima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowAtBorder().allowPlateaus());
        desired[8] = 1.0;
        desired[26] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void extendedLocalMinimum4Test()
    {
        Image res(img);
        res.init(0);

        IImage tmp(img);
	tmp.init(0);

	GridGraph graph(img.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap(static_cast<ViewType>(img)), 
				 out,
				 tmpview,
				 Image::value_type(1.0), // marker
				 NumericTraits<Volume::value_type>::max(), // threshold
				 std::less<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

	//        extendedLocalMinima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);
 
//         res.init(0);
//         localMinima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowPlateaus().neighborhood(4));
//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMinima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().neighborhood(4).allowAtBorder().allowPlateaus());
        desired[8] = 1.0;
        desired[26] = 1.0;
        desired[53] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
   }

    void extendedLocalMaximumTest()
    {
        Image res(img);
        res.init(0);

        IImage tmp(img);
	tmp.init(0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap(static_cast<ViewType>(img)), 
				 out,
				 tmpview,
				 Image::value_type(1.0), // marker
				 NumericTraits<Volume::value_type>::min(), // threshold
				 std::greater<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

//         extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0);

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);
//         res.init(0);
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowPlateaus());
//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowAtBorder().allowPlateaus());
        desired[0] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
   }

    void extendedLocalMaximum4Test()
    {
        Image res(img);
        res.init(0);

        IImage tmp(img);
	tmp.init(0);

	GridGraph graph(img.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap(static_cast<ViewType>(img)), 
				 out,
				 tmpview,
				 Image::value_type(1.0), // marker
				 NumericTraits<Volume::value_type>::min(), // threshold
				 std::greater<Volume::value_type>(),
				 std::equal_to<Volume::value_type>()
				 );

//         extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().allowPlateaus().neighborhood(4));
//         shouldEqualSequence(res.begin(), res.end(), desired);

//         res.init(0);
//         localMaxima(srcImageRange(img), destImage(res), 
//                     LocalMinmaxOptions().neighborhood(4).allowAtBorder().allowPlateaus());
        desired[0] = 1.0;
        desired[27] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);
   }

    void plateauWithHolesTest()
    {
        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 2.0, 3.0, 2.0, 2.0, 1.0, 0.0, 0.0,
            0.0, 1.0, 3.0, 4.0, 4.1, 3.0, 1.0, 0.0, 0.0,
            0.0, 1.0, 2.0, 2.0, 3.0, 2.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        std::copy(in, in+81, img.begin());
        Image res(img.shape(), 0.0);
        IImage tmp(img.shape(), 0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));
	// FIXME: provide more convenient interface!
	extendedLocalMinMaxGraph(graph, 
				 VertexPropMap(static_cast<ViewType>(img)), 
				 out,
				 tmpview,
				 Image::value_type(1.0), // marker
				 NumericTraits<Volume::value_type>::min(), // threshold
				 std::greater<Volume::value_type>(),
				 EqualWithToleranceFunctor<Image::value_type>(0.2)
				 );

//         extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0,
//                             EightNeighborCode(),
//                             EqualWithToleranceFunctor<Image::value_type>(0.2));
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    Image img;
    Volume vol;
};

struct WatershedsTest
{
    typedef MultiArray<2, double> Image;
    typedef MultiArray<2, Int32> IImage;
    typedef Image::difference_type Shape;
    typedef MultiArrayView<2, Image::value_type> ViewType;
    typedef MultiArrayView<2, IImage::value_type> ViewTypeI;

    typedef GridGraphView_CoordsDescriptor<2> GridGraph;
    typedef MultiArrayView_property_map<ViewType> VertexPropMap;
    typedef MultiArrayView_property_map<ViewTypeI> VertexPropMapI;

    WatershedsTest()
	: img(Shape(9,9))
    {
        static const double in[] = {
            0.0,  0.1,  0.1,  0.3,  0.5,  0.3,  0.0,  0.0, 0.0,
            0.0, -0.1,  0.1,  0.0,  1.0,  0.0,  0.3,  0.0, 0.0,
            0.0,  0.5,  2.0,  0.0,  2.0,  2.0,  2.0,  0.0, 0.0,
            0.0,  0.0,  1.0,  1.5,  1.0,  1.0,  3.0,  4.0, 0.0,
            0.3,  0.1,  1.5,  0.0,  0.0,  0.0,  0.0,  3.0, 3.0,
            0.0,  0.0,  0.0,  0.0, -1.0, -1.5, -1.0,  0.0, 0.0,
            0.0,  0.0, -2.0, -2.0, -2.0,  0.0, -2.1, -0.5, 0.0,
            0.0,  0.0,  0.0,  0.0, -1.0,  0.0, -0.1,  0.1, 0.0,
            0.0,  0.0,  0.0,  0.0, -0.5, -0.3, -0.1, -0.1, 0.0};

        ViewType::iterator i = static_cast<ViewType>(img).begin();
        ViewType::iterator end = static_cast<ViewType>(img).end();
        const double * p = in;

        for(; i != end; ++i, ++p)
        {
            // transform data to a range suitable for BucketQueue (in the turbo algorithm)
            *i = *p*10.0 + 30.0;
        }
    }

    void watershedsTest()
    {
        IImage res(img.shape());
        IImage tmp(img.shape(), 0);

	GridGraph graph(img.shape(), IndirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));

        res.init(0);
	
	int count = generateWatershedSeeds_graph(graph,
						 VertexPropMap(static_cast<ViewType>(img)),
						 out,
						 tmpview,
						 SeedOptions().extendedMinima());

        static const double desiredSeeds[] = {
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            0.0,  2.0,  0.0,  3.0,  0.0,  1.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  3.0,  0.0,  0.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  4.0,  4.0,  4.0,  0.0,  5.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};

	shouldEqual(5, count);
        shouldEqualSequence(res.begin(), res.end(), desiredSeeds);

	count = watershedsRegionGrowing_graph(graph,
					      VertexPropMap(static_cast<ViewType>(img)),
					      out,
					      tmpview,
					      WatershedOptions());

        static const double desiredRG[] = {
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  3.0,  1.0,  1.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  1.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0};

	shouldEqual(5, count);
	shouldEqualSequence(res.begin(), res.end(), desiredRG);

        res.init(0);
	count = watershedsRegionGrowing_graph(graph,
					      VertexPropMap(static_cast<ViewType>(img)),
					      out,
					      tmpview,
					      WatershedOptions().keepContours().seedOptions(SeedOptions().extendedMinima()));

        static const double desiredRGC[] = {
            2.0,  2.0,  0.0,  3.0,  0.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  
            0.0,  0.0,  0.0,  4.0,  4.0,  0.0,  5.0,  0.0,  0.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  0.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  0.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0};

        shouldEqual(5, count);
        shouldEqualSequence(res.begin(), res.end(), desiredRGC);

#if 0
	res.init(0);
	count = watershedsRegionGrowing(srcImageRange(img), destImage(res),
					WatershedOptions().turboAlgorithm()
					.seedOptions(SeedOptions().extendedMinima()));
#endif
#if 1
        res.init(0);
	count = watershedsRegionGrowing_graph(graph,
					      VertexPropMap(static_cast<ViewType>(img)),
					      out,
					      tmpview,
					      WatershedOptions().turboAlgorithm()
					      .seedOptions(SeedOptions().extendedMinima()));
#endif

        static const double desiredTRG[] = {
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  5.0,  1.0,  1.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0};

        shouldEqual(5, count);
        shouldEqualSequence(res.begin(), res.end(), desiredTRG);

#if 0
        std::cerr << count << "\n";
        for(int y=0;y<9;++y)
        {
            std::cerr << "            ";
            for(int x=0;x<9;++x)
                std::cerr << res(x,y) << ".0,  ";
            std::cerr << "\n\n";
        }
#endif /* #if 0 */
    }

    void watersheds4Test()
    {
        IImage res(img.shape());
	IImage tmp(img.shape(), 0);

	GridGraph graph(img.shape(), DirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));

        res.init(0);

	int count = generateWatershedSeeds_graph(graph,
						 VertexPropMap(static_cast<ViewType>(img)),
						 out,
						 tmpview,
						 SeedOptions().extendedMinima());

        static const double desiredSeeds[] = {
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            0.0,  2.0,  0.0,  3.0,  0.0,  4.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  3.0,  0.0,  0.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  5.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  6.0,  6.0,  6.0,  0.0,  7.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};

        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredSeeds);

	count = watershedsRegionGrowing_graph(graph,
					      VertexPropMap(static_cast<ViewType>(img)),
					      out,
					      tmpview,
					      WatershedOptions());

        static const double desiredRG[] = {
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  5.0,  5.0,  1.0,  1.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  5.0,  5.0,  5.0,  1.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  5.0,  5.0,  5.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  5.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0};

        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredRG);

        res.init(0);
	count = watershedsRegionGrowing_graph(graph,
					      VertexPropMap(static_cast<ViewType>(img)),
					      out,
					      tmpview,
					      WatershedOptions().keepContours().seedOptions(SeedOptions().extendedMinima()));

        static const double desiredRGC[] = {
            2.0,  2.0,  2.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  4.0,  0.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  0.0,  5.0,  5.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  5.0,  5.0,  5.0,  0.0,  0.0,  
            6.0,  6.0,  6.0,  6.0,  0.0,  5.0,  0.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  0.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  0.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  0.0,  7.0,  7.0};

        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredRGC);

        res.init(0);

#if 1
	count = watershedsRegionGrowing_graph(graph,
						  VertexPropMap(static_cast<ViewType>(img)),
						  out,
						  tmpview,
						  WatershedOptions().turboAlgorithm()
						  .seedOptions(SeedOptions().extendedMinima()));
#endif
#if 0
        count = watershedsRegionGrowing(srcImageRange(img), destImage(res), FourNeighborCode(),
                                        WatershedOptions().turboAlgorithm()
                                          .seedOptions(SeedOptions().extendedMinima()));
#endif

        static const double desiredTRG[] = {
            2.0,  2.0,  2.0,  3.0,  1.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  6.0,  5.0,  7.0,  1.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  5.0,  7.0,  7.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  5.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0};

        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredTRG);

#if 0
        std::cerr << count << "\n";
        for(int y=0;y<9;++y)
        {
            std::cerr << "            ";
            for(int x=0;x<9;++x)
                std::cerr << res(x,y) << ".0,  ";
            std::cerr << "\n";
        }
#endif /* #if 0 */
    }

    Image img;
};

struct RegionGrowingTest
{
    typedef MultiArray<2, double> Image;
    typedef MultiArray<2, Int32> IImage;
    typedef Image::difference_type Shape;
    typedef MultiArrayView<2, Image::value_type> ViewType;
    typedef MultiArrayView<2, IImage::value_type> ViewTypeI;

    typedef GridGraphView_CoordsDescriptor<2> GridGraph;
    typedef MultiArrayView_property_map<ViewType> VertexPropMap;
    typedef MultiArrayView_property_map<ViewTypeI> VertexPropMapI;

    RegionGrowingTest()
	: img(Shape(7,7)), seeds(Shape(7,7))
    {
        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        Image tmp(Shape(7,7));

        ViewType::iterator i = static_cast<ViewType>(tmp).begin();
        ViewType::iterator end = static_cast<ViewType>(tmp).end();
        const double * p = in;

        for(; i != end; ++i, ++p)
        {
            *i = *p;
        }

        distanceTransform(srcImageRange(tmp), destImage(img), 0.0, 2);

        seeds.init(0);

	VertexPropMap seedsView(static_cast<ViewType>(seeds));
	labelGraphWithBackground(GridGraph(tmp.shape(), DirectNeighborhood),
				 VertexPropMap(static_cast<ViewType>(tmp)), 
				 seedsView,
				 0.0);
    }

    struct DirectCostFunctor
    {
        typedef double argument_type;
        typedef double result_type;
        typedef double cost_type;

        void operator()(double const &) {}

        double const & cost(double const & v) const
        {
            return v;
        }
    };

    void voronoiTest()
    {
        Image res(img);

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing(srcImageRange(img), srcImage(seeds),
                            destImage(res), cost);

	ViewType::iterator i = static_cast<ViewType>(res).begin();
        int x,y;

        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = *(i+Shape(x,y));
                double dist1 = VIGRA_CSTD::sqrt((2.0 - x)*(2.0 - x) +
                                         (2.0 - y)*(2.0 - y));
                double dist2 = VIGRA_CSTD::sqrt((5.0 - x)*(5.0 - x) +
                                         (5.0 - y)*(5.0 - y));
                double desired = (dist1 <= dist2) ? 1 : 2;

                if(VIGRA_CSTD::fabs(dist1 - dist2) > 1e-10)
                    shouldEqual(dist, desired);
            }
        }

        IImage wres(img.shape());
	IImage tmp(img.shape(), 0);

	GridGraph graph(img.shape(), DirectNeighborhood);
	// GridGraph graph(img.shape(), IndirectNeighborhood); // same result
	VertexPropMapI out(static_cast<ViewTypeI>(wres));
	VertexPropMapI tmpview(static_cast<ViewTypeI>(tmp));

	watershedsRegionGrowing_graph(graph,
				      VertexPropMap(static_cast<ViewType>(img)),
				      out,
				      tmpview,
				      WatershedOptions().completeGrow()
				      .seedOptions(SeedOptions().minima().threshold(1.0)));

        shouldEqualSequence(res.begin(), res.end(), wres.begin());
    }

    void voronoiWithBorderTest()
    {
        Image res(img);
        Image::value_type reference[] = {
            1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 0,
            1, 1, 1, 1, 1, 0, 2,
            1, 1, 1, 1, 0, 2, 2,
            1, 1, 1, 0, 2, 2, 2,
            1, 1, 0, 2, 2, 2, 2,
            1, 0, 2, 2, 2, 2, 2
        };

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        IImage tmp(img.shape(), 0);

	GridGraph graph(img.shape(), DirectNeighborhood);
	VertexPropMap out(static_cast<ViewType>(res));
	VertexPropMapI tmpRegions(static_cast<ViewTypeI>(tmp));

	seededRegionGrowing_graph(graph,
				  VertexPropMap(static_cast<ViewType>(img)),
				  VertexPropMap(static_cast<ViewType>(seeds)),
				  out,
				  tmpRegions,
				  cost, 
				  KeepContours);

        shouldEqualSequence(res.begin(), res.end(), reference);
    }

    Image img, seeds;
};


struct SimpleAnalysisGridGraphTestSuite
: public vigra::test_suite
{
    SimpleAnalysisGridGraphTestSuite()
    : vigra::test_suite("SimpleAnalysisGridGraphTestSuite")
    {
        add( testCase( &LabelingTest::labelingFourTest1));
        add( testCase( &LabelingTest::labelingFourTest2));
        add( testCase( &LabelingTest::labelingFourTest3));
        add( testCase( &LabelingTest::labelingFourTest4));
        add( testCase( &LabelingTest::labelingToCrackEdgeTest));
        add( testCase( &LabelingTest::labelingEightTest1));
        add( testCase( &LabelingTest::labelingEightTest2));
        add( testCase( &LabelingTest::labelingFourWithBackgroundTest1));
        add( testCase( &LabelingTest::labelingFourWithBackgroundTest2));
        add( testCase( &LabelingTest::labelingEightWithBackgroundTest));

        add( testCase( &LocalMinMaxTest::localMinimumTest));
        add( testCase( &LocalMinMaxTest::localMinimum4Test));
        add( testCase( &LocalMinMaxTest::localMinimumTestThr));
        add( testCase( &LocalMinMaxTest::localMaximumTest));
        add( testCase( &LocalMinMaxTest::localMaximum4Test));
        add( testCase( &LocalMinMaxTest::localMaximumTestThr));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimum4Test));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximum4Test));

        add( testCase( &LocalMinMaxTest::extendedLocalMaximum3DTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimum3DTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximum3DTest2));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimum3DTest2));
        add( testCase( &LocalMinMaxTest::localMaximum3DTest));
        add( testCase( &LocalMinMaxTest::localMinimum3DTest));

        add( testCase( &LocalMinMaxTest::plateauWithHolesTest));
        add( testCase( &WatershedsTest::watershedsTest));
	add( testCase( &WatershedsTest::watersheds4Test));
        add( testCase( &RegionGrowingTest::voronoiTest));
        add( testCase( &RegionGrowingTest::voronoiWithBorderTest));
    }
};

int main(int argc, char ** argv)
{
    SimpleAnalysisGridGraphTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

