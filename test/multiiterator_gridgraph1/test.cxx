#define BGL_COMPAT
#define TEST_BGL_COMPAT

#include <vigra/multi_array.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <vigra/voxelneighborhood.hxx>
#include <vigra/multi_gridgraph.hxx>
#include <vigra/timing.hxx>


using namespace vigra;
using namespace std;

template <class Array, class Shape, class Iter>
void check(Array const & a, Shape const & p, Iter const & i, const char * message)
{
    Shape s = a.shape();
    
    if(p != i.point())
        cerr << message << "coord error at " << i.point() << " (should be " << p << ")\n";
    if(a[p] != (double)i.index())
        cerr << message << "index error at " << p << " (was " << i.index() << ", should be " << a[p] << ")\n";
    AtVolumeBorder atBorder = isAtVolumeBorder(p[0], p[1], p[2], s[0], s[1], s[2]);
    if(atBorder == NotAtBorder)
	{
	    if(i.atBorder())
		cerr << message << "false border at " << p << "\n";
	}
    else
	{
	    if(!i.atBorder())
		cerr << message << "missed border in i at " << p << "\n";
	}
}

int main()
{
    typedef MultiArrayShape<3>::type Shape;
    
    int n = 200;
    Shape s(n, n, n);
    
    MultiArray<3, double> a(s);
    MultiArrayView<3,double> av(a);
    for(int k=0; k<a.size(); ++k)
        a.data()[k] = k;
        
    MultiArray<3, unsigned char> f(s);
    initMultiArrayBorder(destMultiArrayRange(f), 1, 1); // 1 voxel width border around image, is flagged as "1"

    // test three iterators here:
    // one with the flags auxiliary array:
    GridGraphNodeIterator<3, double> i(a.data(), f.data(), a.shape(), a.stride()),
	iend = i.getEndIterator();
    // one without auxiliary array:
    GridGraphNodeIterator2<3, double> i2(a.data(), a.shape(), a.stride()),
	i2end = i2.getEndIterator();
    // and the standard one: a StridedScanOrderIterator:
    MultiArrayView<3,double>::iterator i3 = av.begin(),
	i3end = av.end();


    for(int z=0; z<s[2]; ++z)
	{
	    for(int y=0; y<s[1]; ++y)
		{
		    for(int x=0; x<s[0]; ++x, ++i, ++i2, ++i3)
			{
			    Shape p(x,y,z);
			    // check if each of the iterators get's the current
			    // coordinate right (point==p) as it is scanned over the image:
			    // (also checks the index, which was written as the array data above,
			    // as well as the "i.atBorder()" condition.)
			    check(a, p, i, "i: ");
			    check(a, p, i2, "i2: ");
			    check(a, p, i3, "i3: ");
			}
		}
	}
    cerr << "Tests done\n";
    
    USETICTOC

#define PERFORM_FULL_ARRAY_TESTS
#define USE_PRE_INCREMENTS 1
	// The next tests scan over the image (via ++i oder i+=1) and
	// decrement the voxel value there (after incrementing it in advance by 3
	// to cancel that change if all works well):
#ifdef PERFORM_FULL_ARRAY_TESTS
	TIC
	for(int z=0; z<s[2]; ++z)
	    {
		for(int y=0; y<s[1]; ++y)
		    {
			for(int x=0; x<s[0]; ++x)
			    {
				a(x,y,z) += 3.0;
			    }
		    }
	    }
    TOC
	cout << "<< Direct iteration" << endl;

    i = GridGraphNodeIterator<3, double>(a.data(), f.data(), a.shape(), a.stride());
    TIC
#if USE_PRE_INCREMENTS
	for(; i != iend; ++i)
	    {
		*i -= 1.0;
	    }
#else    
    for(; i != iend; i += 1)
	{
	    *i -= 1.0;
	}
#endif    
    TOC
	cout << "<< GridGraphNodeIterator" << endl;
    i2 = GridGraphNodeIterator2<3, double>(a.data(), a.shape(), a.stride());
    TIC
#if USE_PRE_INCREMENTS
	for(; i2 != i2end; ++i2)
	    {
		*i2 -= 1.0;
	    }
#else    
    for(; i2 != i2end; i2 += 1)
	{
	    *i2 -= 1.0;
	}
#endif    
    TOC
	cout << "<< GridGraphNodeIterator2" << endl;
    i3 = av.begin();
    TIC
#if USE_PRE_INCREMENTS
	for(; i3 != i3end; ++i3)
	    {
		*i3 -= 1.0;
	    }
#else    
    for(; i3 != i3end; i3 += 1)
	{
	    *i3 -= 1.0;
	}
#endif    
    TOC
	cout << "<< original StridedScanOrderIterator" << endl;
#endif



    cout << "Border-Handling iteration test:" << endl;

    // The next tests scan over the image (via ++i oder i+=1) and
    // decrement the voxel value there (after incrementing it in advance by 3
    // to cancel that change if all works well);
    // similar to the above, but sparing the borders.
    TIC
	for(int z=0; z<s[2]; ++z)
	    {
		for(int y=0; y<s[1]; ++y)
		    {
			for(int x=0; x<s[0]; ++x)
			    {
				AtVolumeBorder atBorder = isAtVolumeBorder(x, y, z, s[0], s[1], s[2]);
				if(atBorder == NotAtBorder)
				    a(x,y,z) += 3.0;
			    }
		    }
	    }
    TOC
	cout << "<< Direct iteration" << endl;

    i = GridGraphNodeIterator<3, double>(a.data(), f.data(), a.shape(), a.stride());
    TIC
	for(; i != iend; ++i)
	    {
		bool atBorder = i.atBorder();
		if(!atBorder)
		    *i -= 1.0;
	    }
    TOC
	cout << "<< GridGraphNodeIterator"  << endl;
    
    i2 = GridGraphNodeIterator2<3, double>(a.data(), a.shape(), a.stride());
    TIC
	for(; i2 != i2end; ++i2)
	    {
		bool atBorder = i2.atBorder();
		if(!atBorder)
		    *i2 -= 1.0;
	    }
    TOC
	cout << "<< GridGraphNodeIterator2" << endl;
    
    i3 = av.begin();
    TIC
	for(; i3 != i3end; ++i3)
	    {
		bool atBorder = i3.atBorder();
		if(!atBorder)
		    *i3 -= 1.0;
	    }
    TOC
	cout << "<< original StridedScanOrderIterator" << endl;
    
    // StS: check data again:
    i = GridGraphNodeIterator<3, double>(a.data(), f.data(), a.shape(), a.stride());
    i2 = GridGraphNodeIterator2<3, double>(a.data(), a.shape(), a.stride());
    i3 = av.begin();

    for(int z=0; z<s[2]; ++z)
	{
	    for(int y=0; y<s[1]; ++y)
		{
		    for(int x=0; x<s[0]; ++x, ++i, ++i2, ++i3)
			{
			    Shape p(x,y,z);
			    // check if each of the iterators get's the current
			    // coordinate right (point==p) as it is scanned over the image:
			    // (also checks the index, which was written as the array data above,
			    // as well as the "i.atBorder()" condition.)
			    check(a, p, i, "i: ");
			    check(a, p, i2, "i2: ");
			    check(a, p, i3, "i3: ");
			}
		}
	}
    cerr << "Post-run Tests done\n";
    
    {
	typedef MultiArrayShape<2>::type Shape2;
	ArrayVector<ArrayVector<Shape2> > neighborhood;
	ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
	ArrayVector<ArrayVector<int> > neighborIndexLookup;
	detail::makeArrayNeighborhood(neighborhood, neighborExists, causalNeighborhood, anticausalNeighborhood, neighborIndexLookup, false);

	{
	    ArrayVector<ArrayVector<bool> >::iterator nit = neighborExists.begin(), nend=neighborExists.end();
	    for(;nit!=nend;++nit) {
		ArrayVector<bool>::iterator n2it=nit->begin(), n2end=nit->end();
		for (;n2it!=n2end;++n2it) 
		    cout << (*n2it ? "n":".");
		cout << endl;
	    }
	}

    }

    {
	typedef MultiArrayShape<2>::type Shape2;
	ArrayVector<ArrayVector<Shape2> > neighborhood;
	ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
	ArrayVector<ArrayVector<int> > neighborIndexLookup;

	detail::makeArrayNeighborhood(neighborhood, neighborExists, causalNeighborhood, anticausalNeighborhood, neighborIndexLookup, true);

	{
	    ArrayVector<ArrayVector<bool> >::iterator nit = neighborExists.begin(), nend=neighborExists.end();
	    for(;nit!=nend;++nit) {
		ArrayVector<bool>::iterator n2it=nit->begin(), n2end=nit->end();
		for (;n2it!=n2end;++n2it) 
		    cout << (*n2it ? "n":".");
		cout << endl;
	    }
	}

    }


    { // 3D case?
	typedef MultiArrayShape<3>::type Shape2;
	ArrayVector<ArrayVector<Shape2> > neighborhood;
	ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
	ArrayVector<ArrayVector<int> > neighborIndexLookup;

	detail::makeArrayNeighborhood(neighborhood, neighborExists, causalNeighborhood, anticausalNeighborhood, neighborIndexLookup, false);

	{
	    ArrayVector<ArrayVector<bool> >::iterator nit = neighborExists.begin(), nend=neighborExists.end();
	    for(;nit!=nend;++nit) {
		ArrayVector<bool>::iterator n2it=nit->begin(), n2end=nit->end();
		for (;n2it!=n2end;++n2it) 
		    cout << (*n2it ? "n":".");
		cout << endl;
	    }
	}
	{
	    ArrayVector<ArrayVector<bool> >::iterator nit = causalNeighborhood.begin(), nend=causalNeighborhood.end();
	    for(;nit!=nend;++nit) {
		ArrayVector<bool>::iterator n2it=nit->begin(), n2end=nit->end();
		for (;n2it!=n2end;++n2it) 
		    cout << (*n2it ? "c":".");
		cout << endl;
	    }
	}
	{
	    ArrayVector<ArrayVector<bool> >::iterator nit = anticausalNeighborhood.begin(), nend=anticausalNeighborhood.end();
	    for(;nit!=nend;++nit) {
		ArrayVector<bool>::iterator n2it=nit->begin(), n2end=nit->end();
		for (;n2it!=n2end;++n2it) 
		    cout << (*n2it ? "a":".");
		cout << endl;
	    }
	}

    }

    { // 3D case?
	typedef MultiArrayShape<3>::type Shape2;
	ArrayVector<ArrayVector<Shape2> > neighborhood;
	ArrayVector<ArrayVector<bool> > neighborExists, causalNeighborhood, anticausalNeighborhood;
	ArrayVector<ArrayVector<int> > neighborIndexLookup;

	detail::makeArrayNeighborhood(neighborhood, neighborExists, causalNeighborhood, anticausalNeighborhood, neighborIndexLookup, true);

	{
	    ArrayVector<ArrayVector<bool> >::iterator nit = neighborExists.begin(), nend=neighborExists.end();
	    for(;nit!=nend;++nit) {
		ArrayVector<bool>::iterator n2it=nit->begin(), n2end=nit->end();
		for (;n2it!=n2end;++n2it) 
		    cout << (*n2it ? "n":".");
		cout << endl;
	    }
	}
    }





    { // test new neighbor iterators:
	typedef GridGraphView<3, double> grid_view_type;
	grid_view_type ggv(a);
	grid_view_type::node_iterator it=ggv.get_node_iterator(), itend=ggv.get_node_end_iterator();
	cout << "First node:" << *it << endl;
      
	grid_view_type::neighbor_data_iterator nbit=ggv.get_neighbor_iterator(it), nbend=ggv.get_neighbor_end_iterator(it);
	unsigned int nbnr = 0;
	for (; nbit != nbend; ++nbit,++nbnr) {
	    cout << "neighbor" << nbnr << " has neighborIndex "<< nbit.neighborIndex() << " and dereferences to " << *nbit  << endl;
	    switch (nbnr) {
	    case 0:
		assert(*nbit == 40201);
		assert(nbit.neighborIndex() == 13);
		break;
	    case 1:
		assert(*nbit == 40200);
		assert(nbit.neighborIndex() == 14);
		break;
	    case 2:
		assert(*nbit == 40001);
		assert(nbit.neighborIndex() == 16);
		break;
	    case 3:
		assert(*nbit == 40000);
		assert(nbit.neighborIndex() == 17);
		break;
	    case 4:
		assert(*nbit == 201);
		assert(nbit.neighborIndex() == 22);
		break;
	    case 5:
		assert(*nbit == 200);
		assert(nbit.neighborIndex() == 23);
		break;
	    case 6:
		assert(*nbit == 1);
		assert(nbit.neighborIndex() == 25);
		break;
	    case 7:
		// there's no seventh neighbor
		assert(false);
		break;
	    }
	}
    }


    { // test new neighbor iterators: neighbor_vertex_iterator dereferences to a neighbor (inheriting from StridedScanOrderIterator)
	typedef GridGraphView<3, double> grid_view_type;
	grid_view_type ggv(a);
	grid_view_type::node_iterator it=ggv.get_node_iterator(), itend=ggv.get_node_end_iterator();
	cout << "First node:" << *it << endl;
      
	grid_view_type::neighbor_vertex_iterator nbit=ggv.get_neighbor_vertex_iterator(it), nbend=ggv.get_neighbor_vertex_end_iterator(it);
	unsigned int nbnr = 0;
	for (; nbit != nbend; ++nbit,++nbnr) {
	    cout << "neighbor" << nbnr << ":" << nbit << " has neighborIndex "<< nbit.neighborIndex() << " and dereferences to " << *nbit  << endl;
	    switch (nbnr) {
	    case 0:
		assert(*((grid_view_type::node_iterator)*nbit) == 40201);
		assert(nbit.neighborIndex() == 13);
		break;
	    case 1:
		assert(*((grid_view_type::node_iterator)*nbit) == 40200);
		assert(nbit.neighborIndex() == 14);
		break;
	    case 2:
		assert(*((grid_view_type::node_iterator)*nbit) == 40001);
		assert(nbit.neighborIndex() == 16);
		break;
	    case 3:
		assert(*((grid_view_type::node_iterator)*nbit) == 40000);
		assert(nbit.neighborIndex() == 17);
		break;
	    case 4:
		assert(*((grid_view_type::node_iterator)*nbit) == 201);
		assert(nbit.neighborIndex() == 22);
		break;
	    case 5:
		assert(*((grid_view_type::node_iterator)*nbit) == 200);
		assert(nbit.neighborIndex() == 23);
		break;
	    case 6:
		assert(*((grid_view_type::node_iterator)*nbit) == 1);
		assert(nbit.neighborIndex() == 25);
		break;
	    case 7:
		// there's no seventh neighbor
		assert(false);
		break;
	    }
	}
    }

}
