/************************************************************************/
/*                                                                      */
/*    Copyright 2006-2012 by F. Heinrich, B. Seppke, Stefan Schmidt,    */
/*        Ullrich Koethe                                                */
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
#include <functional>
#include <cmath>
#include "unittest.hxx"

#include <vigra/multi_gridgraph_coords.hxx>
#include <vigra/multi_labelgraph.hxx>

using namespace vigra;
using namespace vigragraph;

struct VolumeLabelingTest
{
    typedef vigra::MultiArray<3,int> IntVolume;
    typedef vigra::MultiArray<3,double> DoubleVolume;

    typedef GridGraphView_CoordsDescriptor<3> GridGraph;
    typedef MultiArrayView<3, IntVolume::value_type> ViewTypeI;
    typedef MultiArrayView<3, DoubleVolume::value_type> ViewTypeD;
    typedef MultiArrayView_property_map<ViewTypeI> VertexPropMapI;
    typedef MultiArrayView_property_map<ViewTypeD> VertexPropMapD;

    VolumeLabelingTest()
    : vol1(IntVolume::difference_type(4,4,4)),vol2(IntVolume::difference_type(4,4,4)),
      vol3(IntVolume::difference_type(5,5,5)),vol4(DoubleVolume::difference_type(5,5,5)),
      vol5(DoubleVolume::difference_type(5,5,5)),vol6(DoubleVolume::difference_type(5,5,5))
    {
        static const int in1[] = { 0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,
                                   0, 0, 0, 0,    0, 1, 1, 0,    0, 1, 1, 0,    0, 0, 0, 0,
                                   0, 0, 0, 0,    0, 1, 1, 0,    0, 1, 1, 0,    0, 0, 0, 0,
                                   0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0,    0, 0, 0, 0};

        IntVolume::iterator i = vol1.begin();
        IntVolume::iterator end = vol1.end();
        const int * p = in1;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

        static const int in2[] = { 0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,
                                   1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1,
                                   0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,
                                   1, 0, 1, 0,    0, 1, 0, 1,    1, 0, 1, 0,    0, 1, 0, 1};

        i = vol2.begin();
        end = vol2.end();
        p = in2;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

                        
        static const int in3[] = { 0, 1, 0, 0, 0,    1, 1, 0, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
                                   1, 1, 1, 0, 0,    1, 1, 1, 0, 0,    1, 1, 1, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
                                   0, 0, 1, 0, 0,    0, 0, 1, 0, 0,    1, 1, 1, 0, 0,    0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
                                   1, 1, 1, 1, 0,    1, 1, 1, 1, 0,    1, 1, 1, 1, 0,    1, 1, 1, 1, 0,    0, 0, 0, 0, 0,
                                   0, 0, 0, 1, 0,    0, 0, 0, 1, 0,    0, 0, 0, 1, 0,    1, 1, 1, 1, 0,    0, 0, 0, 0, 0};

        i = vol3.begin();
        end = vol3.end();
        p = in3;

        for(; i != end; ++i, ++p)
        {
            *i=*p;
        }

        static const double in4[] = { 1.0, 0.0, 0.0, 0.0, 1.0,    0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,    1.0, 0.0, 0.0, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    1.0, 1.0, 1.0, 1.0, 1.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 1.0, 0.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                      1.0, 0.0, 0.0, 0.0, 1.0,    0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 0.0, 1.0, 0.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,    1.0, 0.0, 0.0, 0.0, 1.0};

        DoubleVolume::iterator id = vol4.begin();
        DoubleVolume::iterator endd = vol4.end();
        const double * pd = in4;

        for(; id != endd; ++id, ++pd)
        {
            *id=*pd;
        }

        static const double in5[] = { 0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                      2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 2.0, 2.0, 2.0, 0.0,    0.0, 2.0, 1.0, 2.0, 0.0,    0.0, 2.0, 2.0, 2.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0,
                                      2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,    2.0, 1.0, 0.0, 1.0, 2.0,    2.0, 2.0, 0.0, 2.0, 2.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 1.0, 1.0, 1.0, 0.0,    0.0, 0.0, 0.0, 0.0, 0.0};

        id = vol5.begin();
        endd = vol5.end();
        pd = in5;

        for(; id != endd; ++id, ++pd)
        {
            *id=*pd;
        }

        static const double in6[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 

            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 

            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            1.0, 1.0, 0.0, 1.0, 1.0, 
            1.0, 1.0, 0.0, 1.0, 1.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 

            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            1.0, 1.0, 0.0, 1.0, 1.0, 
            1.0, 1.0, 0.0, 1.0, 1.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 

            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0 };

        id = vol6.begin();
        endd = vol6.end();
        pd = in6;

        for(; id != endd; ++id, ++pd)
        {
            *id=*pd;
        }

    }

    void labelingSixTest1()
    {
        IntVolume res(vol1);
        
	GridGraph graph(vol1.shape(), DirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	unsigned int maxLabel = labelGraph(graph, VertexPropMapI(static_cast<ViewTypeI>(vol1)), out);

        should(2 == maxLabel);

        IntVolume::iterator i1 = vol1.begin();
        IntVolume::iterator i1end = vol1.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == (*i2 - 1.0) );
        }
    }

    void labelingSixTest2()
    {
        IntVolume res(vol2);

	GridGraph graph(vol2.shape(), DirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	should(64 == labelGraph(graph, VertexPropMapI(static_cast<ViewTypeI>(vol2)), out));

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        int address = 0;

        for(; i2 != i2end; ++i2, ++address)
        {
            should( *i2 == address+1 );
        }
    }

    void labelingSixTest3()
    {
        IntVolume res(vol3);

	GridGraph graph(vol3.shape(), DirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
 	should(5 == labelGraph(graph, VertexPropMapI(static_cast<ViewTypeI>(vol3)), out));

        static const int out3[] = { 1, 2, 3, 3, 3,    2, 2, 3, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,
                                    2, 2, 2, 3, 3,    2, 2, 2, 3, 3,    2, 2, 2, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,
                                    4, 4, 2, 3, 3,    4, 4, 2, 3, 3,    2, 2, 2, 3, 3,    3, 3, 3, 3, 3,    3, 3, 3, 3, 3,
                                    2, 2, 2, 2, 3,    2, 2, 2, 2, 3,    2, 2, 2, 2, 3,    2, 2, 2, 2, 3,    3, 3, 3, 3, 3,
                                    5, 5, 5, 2, 3,    5, 5, 5, 2, 3,    5, 5, 5, 2, 3,    2, 2, 2, 2, 3,    3, 3, 3, 3, 3};

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        const int * p = out3;

        for(; i2 != i2end; ++i2, ++p)
        {
            should( *i2 == *p );
        }
    }

    void labelingSixTest4()
    {
        IntVolume res(vol4);

	GridGraph graph(vol4.shape(), DirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
  	should(18 == labelGraph(graph, VertexPropMapD(static_cast<ViewTypeD>(vol4)), out));

        static const int out4[] = { 1, 2, 2, 2, 3,    2, 2, 2, 2, 2,    2, 2, 4, 2, 2,    2, 2, 2, 2, 2,    5, 2, 2, 2, 6,
                                    2, 2, 2, 2, 2,    2, 7, 2, 8, 2,    2, 2, 4, 2, 2,    2, 9, 2,10, 2,    2, 2, 2, 2, 2,
                                    2, 2, 4, 2, 2,    2, 2, 4, 2, 2,    4, 4, 4, 4, 4,    2, 2, 4, 2, 2,    2, 2, 4, 2, 2,
                                    2, 2, 2, 2, 2,    2,11, 2,12, 2,    2, 2, 4, 2, 2,    2,13, 2,14, 2,    2, 2, 2, 2, 2,
                                   15, 2, 2, 2,16,    2, 2, 2, 2, 2,    2, 2, 4, 2, 2,    2, 2, 2, 2, 2,   17, 2, 2, 2,18};

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        const int * p = out4;

        for(; i2 != i2end; ++i2, ++p)
        {
            should( *i2 == *p );
        }
    }

    void labelingSixWithBackgroundTest1()
    {
        IntVolume res(vol5);

	GridGraph graph(vol5.shape(), DirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	unsigned int maxLabel = labelGraphWithBackground(graph, VertexPropMapD(static_cast<ViewTypeD>(vol5)), out, 0);
	should(4 == maxLabel);

        static const int out5[] = { 0, 0, 0, 0, 0,    0, 1, 1, 1, 0,    0, 1, 1, 1, 0,    0, 1, 1, 1, 0,    0, 0, 0, 0, 0,
                                    2, 2, 0, 2, 2,    2, 1, 0, 1, 2,    2, 2, 0, 2, 2,    2, 1, 0, 1, 2,    2, 2, 0, 2, 2,
                                    0, 0, 0, 0, 0,    0, 2, 2, 2, 0,    0, 2, 3, 2, 0,    0, 2, 2, 2, 0,    0, 0, 0, 0, 0,
                                    2, 2, 0, 2, 2,    2, 4, 0, 4, 2,    2, 2, 0, 2, 2,    2, 4, 0, 4, 2,    2, 2, 0, 2, 2,
                                    0, 0, 0, 0, 0,    0, 4, 4, 4, 0,    0, 4, 4, 4, 0,    0, 4, 4, 4, 0,    0, 0, 0, 0, 0};

        IntVolume::iterator i2 = res.begin();
        IntVolume::iterator i2end = res.end();
        const int * p = out5;

        for(; i2 != i2end; ++i2, ++p)
        {
            should( *i2 == *p );
        }
    }

    void labelingTwentySixTest1()
    {
        IntVolume res(vol1);

	GridGraph graph(vol1.shape(), IndirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	should(2 == labelGraph(graph, VertexPropMapI(static_cast<ViewTypeI>(vol1)), out));

        IntVolume::iterator i1 = vol1.begin();
        IntVolume::iterator i1end = vol1.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == (*i2 - 1.0) );
        }
    }

    void labelingTwentySixTest2()
    {
        IntVolume res(vol2);

	GridGraph graph(vol2.shape(), IndirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	should(2 == labelGraph(graph, VertexPropMapI(static_cast<ViewTypeI>(vol2)), out));

        IntVolume::iterator i1 = vol2.begin();
        IntVolume::iterator i1end = vol2.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == (*i2 - 1.0) );
        }
    }

    void labelingTwentySixTest3()
    {
        IntVolume res(vol4);

	GridGraph graph(vol4.shape(), IndirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
	should(2 == labelGraph(graph, VertexPropMapD(static_cast<ViewTypeD>(vol4)), out));

        DoubleVolume::iterator i1 = vol4.begin();
        DoubleVolume::iterator i1end = vol4.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == 2-*i2 );
        }
    }

    void labelingTwentySixWithBackgroundTest1()
    {
        IntVolume res(vol5);

	GridGraph graph(vol5.shape(), IndirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
 	should(2 == labelGraphWithBackground(graph, VertexPropMapD(static_cast<ViewTypeD>(vol5)), out, 0));

        DoubleVolume::iterator i1 = vol5.begin();
        DoubleVolume::iterator i1end = vol5.end();
        IntVolume::iterator i2 = res.begin();

        for(; i1 != i1end; ++i1, ++i2)
        {
            should( *i1 == *i2 );
        }
    }

    void labelingAllTest()
    {
        IntVolume res(vol6.shape());
        static const int out6[] = {
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 

                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 

                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                1, 1, 0, 2, 2, 
                1, 1, 0, 2, 2, 
                0, 0, 0, 0, 0, 

                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                1, 1, 0, 2, 2, 
                1, 1, 0, 2, 2, 
                0, 0, 0, 0, 0, 

                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0, 
                0, 0, 0, 0, 0 };


	GridGraph graph(vol6.shape(), IndirectNeighborhood);
	GridGraph graphDirect(vol6.shape(), DirectNeighborhood);
	VertexPropMapI out(static_cast<ViewTypeI>(res));
 	should(2 == labelGraphWithBackground(graphDirect, VertexPropMapD(static_cast<ViewTypeD>(vol6)), out, 0));

        shouldEqualSequence(res.begin(), res.end(), out6);

 	should(2 == labelGraphWithBackground(graph, VertexPropMapD(static_cast<ViewTypeD>(vol6)), out, 0));
        shouldEqualSequence(res.begin(), res.end(), out6);

 	should(3 == labelGraph(graphDirect, VertexPropMapD(static_cast<ViewTypeD>(vol6)), out));
        res -= 1;
        shouldEqualSequence(res.begin(), res.end(), out6);

 	should(3 == labelGraph(graph, VertexPropMapD(static_cast<ViewTypeD>(vol6)), out));
        res -= 1;
        shouldEqualSequence(res.begin(), res.end(), out6);

    }

    IntVolume vol1, vol2, vol3;
    DoubleVolume vol4, vol5, vol6;
};



struct VolumeLabelingTestSuite
: public vigra::test_suite
{
    VolumeLabelingTestSuite()
    : vigra::test_suite("VolumeLabelingTestSuite")
    {
        add( testCase( &VolumeLabelingTest::labelingSixTest1));
        add( testCase( &VolumeLabelingTest::labelingSixTest2));
        add( testCase( &VolumeLabelingTest::labelingSixTest3));
        add( testCase( &VolumeLabelingTest::labelingSixTest4));
        add( testCase( &VolumeLabelingTest::labelingSixWithBackgroundTest1));
        add( testCase( &VolumeLabelingTest::labelingTwentySixTest1));
        add( testCase( &VolumeLabelingTest::labelingTwentySixTest2));
        add( testCase( &VolumeLabelingTest::labelingTwentySixTest3));
        add( testCase( &VolumeLabelingTest::labelingTwentySixWithBackgroundTest1));
        add( testCase( &VolumeLabelingTest::labelingAllTest));
    }
};

int main(int argc, char ** argv)
{
    VolumeLabelingTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}

