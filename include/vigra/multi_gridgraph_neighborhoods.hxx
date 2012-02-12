#ifndef VIGRA_MULTI_GRIDGRAPH_NEIGHBORHOODS_HXX
#define VIGRA_MULTI_GRIDGRAPH_NEIGHBORHOODS_HXX

#include "array_vector.hxx"

namespace vigra {

    enum NeighborhoodType { IndirectNeighborhood=false, DirectNeighborhood=true };

namespace detail {

template <unsigned int Level>
struct MakeDirectArrayNeighborhood
{
    template <class Array>
    static void points(Array & a)
    {
        typedef typename Array::value_type Shape;
        
        Shape point;
        point[Level] = -1;
        a.push_back(point);
        MakeDirectArrayNeighborhood<Level-1>::points(a);
        point[Level] = 1;
        a.push_back(point);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int i)
    {
        a.push_back((i & (1 << 2*Level)) == 0);
        MakeDirectArrayNeighborhood<Level-1>::exists(a, i);
        a.push_back((i & (2 << 2*Level)) == 0);
    }
};

template <>
struct MakeDirectArrayNeighborhood<0>
{
    template <class Array>
    static void points(Array & a)
    {
        typedef typename Array::value_type Shape;
        
        Shape point;
        point[0] = -1;
        a.push_back(point);
        point[0] = 1;
        a.push_back(point);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int i)
    {
        a.push_back((i & 1) == 0);
        a.push_back((i & 2) == 0);
    }
};

template <unsigned int Level>
struct MakeIndirectArrayNeighborhood
{
    template <class Array, class Shape>
    static void points(Array & a, Shape point, bool isCenter = true)
    {
        point[Level] = -1;
        MakeIndirectArrayNeighborhood<Level-1>::points(a, point, false);
        point[Level] = 0;
        MakeIndirectArrayNeighborhood<Level-1>::points(a, point, isCenter);
        point[Level] = 1;
        MakeIndirectArrayNeighborhood<Level-1>::points(a, point, false);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int i, bool isCenter = true)
    {
        if((i & (1 << 2*Level)) == 0)
        {
            MakeIndirectArrayNeighborhood<Level-1>::exists(a, i, false);
        } else {
	  // case was missing!! 
	  // need to recurse to flag all those downstream NB as missing, 
	  // or push the corresponding number of "falses"
	  MakeIndirectArrayNeighborhood<Level-1>::existsNegativeEnumerate(a);
	}
        MakeIndirectArrayNeighborhood<Level-1>::exists(a, i, isCenter);
        if((i & (2 << 2*Level)) == 0)
        {
            MakeIndirectArrayNeighborhood<Level-1>::exists(a, i, false);
        } else {
	  // case was missing!! 
	  // need to recurse to flag all those downstream NB as missing, 
	  // or push the corresponding number of "falses"
	  MakeIndirectArrayNeighborhood<Level-1>::existsNegativeEnumerate(a);
	}
    }

    template <class Array>
    static void existsNegativeEnumerate(Array & a)
    {
      MakeIndirectArrayNeighborhood<Level-1>::existsNegativeEnumerate(a);
      MakeIndirectArrayNeighborhood<Level-1>::existsNegativeEnumerate(a);
      MakeIndirectArrayNeighborhood<Level-1>::existsNegativeEnumerate(a);
    }

};

template <>
struct MakeIndirectArrayNeighborhood<0>
{
    template <class Array, class Shape>
    static void points(Array & a, Shape point, bool isCenter = true)
    {
        point[0] = -1;
        a.push_back(point);
        if(!isCenter) // the center point is not a neighbor, it's just convenient to do the enumeration this way...
        {
            point[0] = 0;
            a.push_back(point);
        }
        point[0] = 1;
        a.push_back(point);
    }
    
    template <class Array>
    static void exists(Array & a, unsigned int i, bool isCenter = true)
    {
        a.push_back((i & 1) == 0);
        if(!isCenter)
        {
            a.push_back(true);
        }
        a.push_back((i & 2) == 0);
    }

    template <class Array>
    static void existsNegativeEnumerate(Array & a)
    {
      a.push_back(false);
      a.push_back(false);
      a.push_back(false);
    }

};

template <class Array>
void orderNeighborhood2D(Array & a, 
                         bool directNeighborsOnly = true)
{
    typename Array::value_type v = a[0];
    if(directNeighborsOnly)
    {
        a[0] = a[2];
        a[2] = a[1];
        a[1] = v;
    }
    else
    {
        a[0] = a[4];
        a[4] = a[3];
        a[3] = v;
        v = a[2];
        a[2] = a[1];
        a[1] = v;
    }
}

template <class Shape>
void
makeArrayNeighborhood(ArrayVector<ArrayVector<Shape> > & neighborOffsets, 
                      ArrayVector<ArrayVector<bool> > & neighborExists,
                      ArrayVector<ArrayVector<bool> > & causalNeighborExists,
                      ArrayVector<ArrayVector<bool> > & anticausalNeighborExists,
		      ArrayVector<ArrayVector<int> > & neighborIndexLookup,
                      bool directNeighborsOnly = true)
{
    enum { N = Shape::static_size };
    unsigned int size = 1 << 2*N;
    Shape strides = cumprod(Shape(MultiArrayIndex(3))) / 3; 
    
    neighborOffsets.resize(size);
    neighborOffsets[0].clear(); // [0] is the standard case of all neighbors present
    if(directNeighborsOnly)
    {
        MakeDirectArrayNeighborhood<N-1>::points(neighborOffsets[0]);
    }
    else
    {
      Shape point; // represents the center
        MakeIndirectArrayNeighborhood<N-1>::points(neighborOffsets[0], point);
    }
    
    unsigned int neighborCount = neighborOffsets[0].size(); // maximal number of neighbors

#ifdef VERBOSE    
    std::cerr << " size " << neighborCount << ": " << neighborOffsets[0] << "\n strides ";
    for(unsigned int l=0; l<neighborCount; ++l)
        std::cerr << dot(neighborOffsets[0][l], strides) << ", ";
    std::cerr << "\n\n";
#endif
    
    neighborExists.resize(size);
    causalNeighborExists.resize(size);
    anticausalNeighborExists.resize(size);
    neighborIndexLookup.resize(size);

    for(unsigned int k=0; k<size; ++k) // iterate all k neighborhood codes
    {
	if (k>0) 
	    neighborOffsets[k].clear();
        neighborExists[k].clear();
        if(directNeighborsOnly)
        {
            MakeDirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
        else
        {
            MakeIndirectArrayNeighborhood<N-1>::exists(neighborExists[k], k);
        }
        
        causalNeighborExists[k].resize(neighborCount);
        anticausalNeighborExists[k].resize(neighborCount);
        
        for(unsigned int l = 0; l<neighborCount; ++l)
        {
            MultiArrayIndex stride = dot(neighborOffsets[0][l], strides);
            if(stride < 0)
            {
                causalNeighborExists[k][l] = neighborExists[k][l];
                anticausalNeighborExists[k][l] = false;
            }
            else
            {
                causalNeighborExists[k][l] = false;
                anticausalNeighborExists[k][l] = neighborExists[k][l];
            }
	    if (neighborExists[k][l])
		neighborIndexLookup[k].push_back(l);
	    if (k>0)
		if (neighborExists[k][l])
		    neighborOffsets[k].push_back(neighborOffsets[0][l]);
        }
    }

}




template <class Shape>
void
makeArraySubNeighborhood(const ArrayVector<Shape> & allNeighborOffsets, 
			 const ArrayVector<ArrayVector<bool> > & neighborExists,
			 const Shape strides,
			 ArrayVector<ArrayVector<MultiArrayIndex> > & neighborIndices
			 )
{
    enum { N = Shape::static_size };
    unsigned int size = 1 << 2*N;
    
    neighborIndices.resize(size);
    const unsigned int neighborCount = allNeighborOffsets.size(); // maximal number of neighbors
    
    for (unsigned int k=0; k<size; ++k)  // iterate all k neighborhood codes
	for(unsigned int l=0; l<neighborCount; ++l) 
	    if (neighborExists[k][l])
		neighborIndices[k].push_back(dot(allNeighborOffsets[l], strides));
#if 0
    for (unsigned int k=0; k<size; ++k)  // iterate all k neighborhood codes
    {
	std::cerr << " NB-type " << k << ": ";
	for(unsigned int l=0; l<neighborCount; ++l) 
	    if (neighborExists[k][l])
	    {
		std::cerr << neighborIndices[k].back() << ", ";
	    }
	std::cerr << std::endl;
    }
#endif
}


} // namespace vigra::detail

} // namespace vigra


#endif // VIGRA_MULTI_GRIDGRAPH_NEIGHBORHOODS_HXX
