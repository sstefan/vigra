#include <typeinfo>
#include <iostream>

#include <unittest.hxx>


#include <vigra/timing.hxx>

#include <vigra/multi_array.hxx>
#include <vigra/multi_iterator_coupled.hxx>

#include "multi_testdata.hxx"


USETICTOC

using namespace std;
using namespace vigra;

const bool verbose = true;
const size_t repetitions = 1;
const size_t best_of_repetitions = 3;



struct CoupledMultiImageIteratorTests {
    enum dimension { dim=3 };
#if 1
    typedef TestData3_nD<unsigned char, dim> tTestData1;
    typedef TestData3_nD<unsigned int, dim> tTestData2;
    typedef TestData3_nD<double, dim> tTestData3;
#else
    typedef TestData3_nD<double, dim> tTestData1;
    typedef TestData3_nD<double, dim> tTestData2;
    typedef TestData3_nD<double, dim> tTestData3;
#endif
    typedef tTestData1::view_type tv1;
    typedef tTestData2::view_type tv2;
    typedef tTestData3::view_type tv3;
    
    tTestData1 td1;
    tTestData2 td2;
    tTestData3 td3;
    
    tv1 &i1;
    tv2 &i2;
    tv3 &i3;

    CoupledMultiImageIteratorTests()
    : i1(td1.vol_view),
      i2(td2.vol_view),
      i3(td3.vol_view) 
    {}



    void coupledIteratorTestChained()
    {
	
	typedef ChainedPtrHolder<tv1, 
	    ChainedPtrHolder<tv3>::type >::type Chain;

	typedef CoupledScanOrderIterator<dim, Chain> 
					 coupled_iter_type;
	
 	coupled_iter_type it(i1.shape(), 
 			     Chain(i1.data(), i1.stride(), 
 				   Chain::next_type(i3.data(), i3.stride()))),
 	    itend = it.getEndIterator();

 	std::cerr << "Coupled Iterator size=" << sizeof(it)  << std::endl;
	double cnt=0;
	for (; it != itend; ++it) {
	    *it.top().next() = *it.top() * *it.top(); // square first image, store in third
	    cnt += *(it.top());
	}
	shouldEqual(cnt, 630070);

	// check with separate iteration
	tTestData1::view_type::iterator i1 = td1.vol_view.begin(), ie1=td1.vol_view.end();
	tTestData3::view_type::iterator i3 = td3.vol_view.begin();
	for (;i1 != ie1; ++i1, ++i3) {
	    shouldEqual(*i3, *i1 * *i1);
	}
    }

    //! another version with a nicer interface
    void coupledIteratorTestChainedNicerInterface()
    {
	
	typedef ChainedPtrHolder<tv1, 
	    ChainedPtrHolder<tv3>::type >::type Chain;

	typedef CoupledScanOrderIterator<dim, Chain> 
					 coupled_iter_type;
	
	coupled_iter_type it(i1.shape(), 
			     Chain(i1.data(), i1.stride(), 
				   Chain::next_type(i3.data(), i3.stride()))),
	    itend = it.getEndIterator();

	for (; it != itend; ++it) {
	    *it.get<1>() = *it.get<0>() * *it.get<0>(); // square first image, store in third	    
	}

	// check with separate iteration
	tTestData1::view_type::iterator i1 = td1.vol_view.begin(), ie1=td1.vol_view.end();
	tTestData3::view_type::iterator i3 = td3.vol_view.begin();
	for (;i1 != ie1; ++i1, ++i3) {
	    shouldEqual(*i3, *i1 * *i1);
	}
    }



    void coupledIteratorTestFactory_AddImages()
    {
	typedef CoupledScanOrderIteratorFactory<dim, tv1, tv2, tv3> IterFactory;
	typedef IterFactory::coupled_iter_type coupled_iter_type;
	
	coupled_iter_type it = IterFactory::makeCoupledIterator(i1, i2, i3),
	    itend = it.getEndIterator();
	
	for (; it != itend; ++it) {
	    // what kind of interface to use here?
	    //  it[] is already taken; operator() ?
	    //			*it(2) = *it(0)+*it(1);
	    *it.get<2>() = *it.get<0>() + *it.get<1>();
	    //	    *it.component3() = *it.component1() + *it.component2();
	}

	// check with separate iteration
	tTestData1::view_type::iterator i1 = td1.vol_view.begin(), ie1=td1.vol_view.end();
	tTestData2::view_type::iterator i2 = td2.vol_view.begin();
	tTestData3::view_type::iterator i3 = td3.vol_view.begin();
	for (;i1 != ie1; ++i1,++i2, ++i3) {
	    shouldEqual(*i3, *i1+*i2);
	}
    }



    void coupledIteratorTestFactory_ImagePair()
    {
	typedef CoupledScanOrderIteratorFactory<dim, tv1, tv3> IterFactory;
	typedef IterFactory::coupled_iter_type coupled_iter_type;
	
	coupled_iter_type it = IterFactory::makeCoupledIterator(i1, i3),
	    itend = it.getEndIterator();
	
	for (; it != itend; ++it) {
	    *it.get<1>() = *it.get<0>() * *it.get<0>(); // square first image, store in third
	}

	// check with separate iteration
	tTestData1::view_type::iterator i1 = td1.vol_view.begin(), ie1=td1.vol_view.end();
	tTestData3::view_type::iterator i3 = td3.vol_view.begin();
	for (;i1 != ie1; ++i1, ++i3) {
	    shouldEqual(*i3, *i1 * *i1);
	}
    }


    //! Test the coupled iterator with 0 associated images, which just yields the coordinates.
    void coupledIteratorTestFactory_Nullary()
    {
	typedef CoupledScanOrderIteratorFactory<dim> IterFactory;
	typedef IterFactory::coupled_iter_type coupled_iter_type;
	
	coupled_iter_type it(i1.shape()),
	    itend = it.getEndIterator();

	for (; it != itend; ++it) {
	    shouldEqual( (*it)[0], it.point()[0] );
	    shouldEqual( (*it)[1], it.point()[1] );
	    shouldEqual( (*it)[2], it.point()[2] );
	}
    }


    //! Test the coupled iterator with 0 associated images, which just yields the coordinates.
    void coupledIteratorTestFactory_Nullary_ComponentAccess()
    {
	typedef CoupledScanOrderIteratorFactory<dim> IterFactory;
	typedef IterFactory::coupled_iter_type coupled_iter_type;
	
	coupled_iter_type it(i1.shape()),
	    itend = it.getEndIterator();

	for (; it != itend; ++it) {
	    shouldEqual( (*it.get<0>())[0], it.point()[0] );
	    shouldEqual( (*it.get<0>())[1], it.point()[1] );
	    shouldEqual( (*it.get<0>())[2], it.point()[2] );
	}
    }



    //! "virtual components" (those with index >= number of
    //  actual bands) are mapped to the iterator itself
    //  and derefence as the coordinate tuple as well
    void coupledIteratorTestFactory_ImagePair_CheckVirtualComponents()
    {
	
	typedef ChainedPtrHolder<tv1, 
	    ChainedPtrHolder<tv3>::type >::type Chain;

	typedef CoupledScanOrderIterator<dim, Chain> 
					 coupled_iter_type;
	
	coupled_iter_type it(i1.shape(), 
			     Chain(i1.data(), i1.stride(), 
				   Chain::next_type(i3.data(), i3.stride()))),
	    itend = it.getEndIterator();

	tTestData1::view_type::iterator i1 = td1.vol_view.begin(), ie1=td1.vol_view.end();
	tTestData3::view_type::iterator i3 = td3.vol_view.begin();

	for (; it != itend; ++it, ++i1, ++i3) {
	    shouldEqual( &(*it.get<0>()), &(*i1) );
	    shouldEqual( &(*it.get<1>()), &(*i3) );

	    shouldEqual( (*it)[0], it.point()[0] );
	    shouldEqual( (*it)[1], it.point()[1] );
	    shouldEqual( (*it)[2], it.point()[2] );

	    shouldEqual( (*it.get<2>())[0], it.point()[0] );
	    shouldEqual( (*it.get<2>())[1], it.point()[1] );
	    shouldEqual( (*it.get<2>())[2], it.point()[2] );

	    shouldEqual( (*it.get<3>())[0], it.point()[0] );
	    shouldEqual( (*it.get<3>())[1], it.point()[1] );
	    shouldEqual( (*it.get<3>())[2], it.point()[2] );

	    shouldEqual( (*it.get<15>())[0], it.point()[0] );
	    shouldEqual( (*it.get<15>())[1], it.point()[1] );
	    shouldEqual( (*it.get<15>())[2], it.point()[2] );
	}
    }



    void coupledIteratorTestFactory_Unary()
    {
	typedef CoupledScanOrderIteratorFactory<dim, tv3> IterFactory;
	typedef IterFactory::coupled_iter_type coupled_iter_type;

	coupled_iter_type it = IterFactory::makeCoupledIterator(i3),
	    itend = it.getEndIterator();
	
	for (; it != itend; ++it) {
	    *it.get<0>() = it.point()[0] + it.shape()[0]*it.point()[1] + it.shape()[0]*it.shape()[1]*it.point()[2];
	}

	// check with separate iteration
	tTestData3::view_type::iterator i3 = td3.vol_view.begin(), ie3=td3.vol_view.end();
	for (;i3 != ie3; ++i3) {
	    shouldEqual(*i3, i3.point()[0] + i3.shape()[0]*i3.point()[1] + i3.shape()[0]*i3.shape()[1]*i3.point()[2]);
	}
    }



    //! Test coupled iterator creation by factory 
    void coupledIteratorTestChained2()
    {
	typedef CoupledScanOrderIteratorFactory<dim, tv1, tv3> IterFactory;

	typedef IterFactory::coupled_iter_type coupled_iter_type;
	
	coupled_iter_type it = IterFactory::makeCoupledIterator(i1, i3),
	    itend = it.getEndIterator();

	for (; it != itend; ++it) {
	    *it.top().next() = *it.top() * *it.top(); // square first image, store in third
	}

	// check with separate iteration
	tTestData1::view_type::iterator i1 = td1.vol_view.begin(), ie1=td1.vol_view.end();
	tTestData3::view_type::iterator i3 = td3.vol_view.begin();
	for (;i1 != ie1; ++i1, ++i3) {
	    shouldEqual(*i3, *i1 * *i1);
	}
    }



};



struct CoupledMultiImageIteratorTests2 {
    typedef MultiArray <1, unsigned char> array1_t;
    typedef array1_t::difference_type shape1_t;
    typedef MultiArray <3, unsigned char> array3_t;
    typedef array3_t::difference_type shape3_t;
    typedef MultiArrayView<3, unsigned char> array_view3_t;
    typedef array_view3_t::iterator iterator3_t;
    
    typedef CoupledScanOrderIteratorFactory<3, array_view3_t> IterFactory;
    typedef IterFactory::coupled_iter_type coupled_iter_type;	

    shape3_t s;
    array3_t a3;
    
    CoupledMultiImageIteratorTests2()
    : s(shape3_t(2,3,5)),
      a3(s, 1)
    {
    }


    void test_iterator ()
    {
        // test scan-order navigation
        array_view3_t av = a3;
	coupled_iter_type
	    i1 = IterFactory::makeCoupledIterator(av),
	    i2 = IterFactory::makeCoupledIterator(av),
	    i3 = IterFactory::makeCoupledIterator(av),
	    iend = i1.getEndIterator();

#if 0
	// these are not directly possible for coupled iterator:
        shouldEqual(&i1[0], &a3(0,0,0));
        shouldEqual(&i1[1], &a3(1,0,0));
        shouldEqual(&i1[2], &a3(0,1,0));
        shouldEqual(&i1[3], &a3(1,1,0));
        shouldEqual(&i1[6], &a3(0,0,1));
        shouldEqual(&i1[7], &a3(1,0,1));
        shouldEqual(&i1[9], &a3(1,1,1));
#endif

        shouldEqual(&*(i1+0).top(), &a3(0,0,0));
        shouldEqual(&*(i1+1).top(), &a3(1,0,0));
        shouldEqual(&*(i1+2).top(), &a3(0,1,0));
        shouldEqual(&*(i1+3).top(), &a3(1,1,0));
        shouldEqual(&*(i1+6).top(), &a3(0,0,1));
        shouldEqual(&*(i1+7).top(), &a3(1,0,1));
        shouldEqual(&*(i1+9).top(), &a3(1,1,1));

        shouldEqual(&*(i1+shape3_t(0,0,0)).top(), &a3(0,0,0));
        shouldEqual(&*(i1+shape3_t(1,0,0)).top(), &a3(1,0,0));
        shouldEqual(&*(i1+shape3_t(0,1,0)).top(), &a3(0,1,0));
        shouldEqual(&*(i1+shape3_t(1,1,0)).top(), &a3(1,1,0));
        shouldEqual(&*(i1+shape3_t(0,0,1)).top(), &a3(0,0,1));
        shouldEqual(&*(i1+shape3_t(1,0,1)).top(), &a3(1,0,1));
        shouldEqual(&*(i1+shape3_t(1,1,1)).top(), &a3(1,1,1));

        shouldEqual(&*(iend-1).top(), &a3(1,2,4));
        shouldEqual(&*(iend-2).top(), &a3(0,2,4));
        shouldEqual(&*(iend-3).top(), &a3(1,1,4));
        shouldEqual(&*(iend-7).top(), &a3(1,2,3));
        shouldEqual(&*(iend-8).top(), &a3(0,2,3));
        shouldEqual(&*(iend-10).top(), &a3(0,1,3));

	i3 = iend-1;
	shouldEqual(&*(i3-shape3_t(0,0,0)).top(), &a3(1,2,4));
	shouldEqual(&*(i3-shape3_t(1,0,0)).top(), &a3(0,2,4));
        shouldEqual(&*(i3-shape3_t(0,1,0)).top(), &a3(1,1,4));
        shouldEqual(&*(i3-shape3_t(1,1,0)).top(), &a3(0,1,4));
	shouldEqual(&*(i3-shape3_t(0,0,1)).top(), &a3(1,2,3));
	shouldEqual(&*(i3-shape3_t(1,0,1)).top(), &a3(0,2,3));
	shouldEqual(&*(i3-shape3_t(1,1,1)).top(), &a3(0,1,3));

#if 0
	// These are not directly possible
	// for the coupled iterator:
        shouldEqual(&iend[-1], &a3(1,2,4));
        shouldEqual(&iend[-2], &a3(0,2,4));
        shouldEqual(&iend[-3], &a3(1,1,4));
        shouldEqual(&iend[-7], &a3(1,2,3));
        shouldEqual(&iend[-8], &a3(0,2,3));
        shouldEqual(&iend[-10], &a3(0,1,3));
#endif

	i3 = i1;
	i3 += shape3_t(0,0,1);
	shouldEqual(i3.index(), 6);
	shouldEqual(i3.point(), shape3_t(0,0,1));
	i3 -= shape3_t(0,0,1);
	shouldEqual(i3.index(), 0);
	shouldEqual(i3.point(), shape3_t(0,0,0));
	should(i3 == i1);

        unsigned int count = 0;
        shape3_t p;

        // iterate over the third dimension
        for (p[2]=0; p[2] != s[2]; ++p[2]) 
        {
            for (p[1]=0; p[1] != s[1]; ++p[1]) 
            {
                for (p[0]=0; p[0] != s[0]; ++p[0], ++i1, i2 += 1, ++count)
                {
                    shouldEqual(&*i1.top(), &a3[p]);
                    shouldEqual(&*i2.top(), &a3[p]);
                    shouldEqual(i1.top().operator->(), &a3[p]);
                    shouldEqual(i2.top().operator->(), &a3[p]);
                    shouldEqual(i1.point(), p);
                    shouldEqual(i2.point(), p);
                    shouldEqual(i1.index(), count);
                    shouldEqual(i2.index(), count);

                    should(i1 != iend);
                    should(!(i1 == iend));
                    should(i1 < iend);
                    should(i1 <= iend);
                    should(!(i1 > iend));
                    should(!(i1 >= iend));

                    shouldEqual(iend - i1, av.size() - count);

                    bool atBorder = p[0] == 0 || p[0] == s[0]-1 || p[1] == 0 || p[1] == s[1]-1 ||
                                    p[2] == 0 || p[2] == s[2]-1;
                    if(!atBorder)
                    {
                        should(!i1.atBorder());
                        should(!i2.atBorder());
                    }
                    else
                    {
                        should(i1.atBorder());
                        should(i2.atBorder());
                    }
                }
            }
        }

        should(i1 == iend);
        should(!(i1 != iend));
        should(!(i1 < iend));
        should(i1 <= iend);
        should(!(i1 > iend));
        should(i1 >= iend);

        should(i2 == iend);
        should(!(i2 != iend));
        should(!(i2 < iend));
        should(i2 <= iend);
        should(!(i2 > iend));
        should(i2 >= iend);

        shouldEqual(iend - i1, 0);
        shouldEqual(iend - i2, 0);
        shouldEqual (count, av.size());

	--i1;
        i2 -= 1;
        shouldEqual(&*i1.top(), &a3(1,2,4));
        shouldEqual(&*i2.top(), &a3(1,2,4));
    }


};


struct GridGraphsTestSuite
: public vigra::test_suite
{
    GridGraphsTestSuite()
    : vigra::test_suite("GridGraphsTestSuite")
    {
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestChained));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestChainedNicerInterface));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestFactory_AddImages));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestFactory_ImagePair));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestFactory_Nullary));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestFactory_Nullary_ComponentAccess));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestFactory_ImagePair_CheckVirtualComponents));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestFactory_Unary));
	add( testCase( &CoupledMultiImageIteratorTests::coupledIteratorTestChained2));

	add( testCase( &CoupledMultiImageIteratorTests2::test_iterator));
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
