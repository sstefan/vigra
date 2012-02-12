// Testing code: Adapting to grid-graphs and general graphs
//#define TEST_ONLY

/************************************************************************/
/*                                                                      */
/*    Copyright 1998-2012 by Stefan Schmidt, Ullrich Koethe             */
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


#ifndef VIGRA_MULTI_LOCALMINMAX_HXX
#define VIGRA_MULTI_LOCALMINMAX_HXX

#include <vector>
#include <functional>
#include <vigra/multi_array.hxx> 
#include <vigra/graphs.hxx>
#include <vigra/localminmax.hxx>
#include <vigra/multi_labelgraph.hxx>


namespace vigra {


    //! localMinMax() variant for BGL-like interfaces
    template <class Graph, class T1Map, class T2Map, class Compare>
void
localMinMaxGraph(
	     Graph &G, 
 	     T1Map &src,
 	     T2Map &dest,
	     typename vigragraph::property_traits<T2Map>::value_type marker, // unsigned int neighborhood,
	     typename vigragraph::property_traits<T1Map>::value_type threshold,
	     Compare compare)
{
    typedef typename vigragraph::graph_traits<Graph>::vertex_iterator graph_scanner;
    typedef typename vigragraph::graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename vigragraph::property_traits<T1Map>::value_type T1;
    typedef typename vigragraph::property_traits<T2Map>::value_type T2;

    graph_scanner srcit, srcend;
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(G);

    neighbor_iterator nbit, nbend;
#ifdef TEST_ONLY    
    size_t counter=0;
    size_t counter1=0;
    T1 dummysum=0;
#endif
    for (; srcit != srcend; ++srcit) {
#ifdef TEST_ONLY    
      ++counter1;
#endif
	  const T1 refval = src[*srcit];

	  if (!compare(refval, threshold))
	      continue;
	  
	  vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*srcit, G);
	  
	  bool local_extremum = true;
	  int j=0;
	  for (;nbit != nbend; ++nbit) {
	    ++j;
#ifdef TEST_ONLY    
	    ++counter;
#endif
	    if (!compare(refval, src[*nbit])) {
	      local_extremum = false;
	      break;
	    }
	  }
	  if (local_extremum)
#ifdef TEST_ONLY    
	      dummysum += marker;
#else	 
  	      dest[*srcit] = marker;
#endif

	  // TODO: case of no neighbors => local min/max per definition?
    }
#ifdef TEST_ONLY    
    std::cout << " counter=" << counter << " counter1=" << counter1 << std::endl;
    std::cout << " dummysum=" << ((int)dummysum) << std::endl;
#endif

}





  // Two-Iterator-Variant of the above:
    template <class Graph, class Graph2, class T1Map, class T2Map, class Compare>
void
localMinMaxGraph2(
	     Graph &Gsrc, 
 	     T1Map &src,
	     Graph2 &Gdest, 
 	     T2Map &dest,
	     typename vigragraph::property_traits<T2Map>::value_type marker, // unsigned int neighborhood,
	     typename vigragraph::property_traits<T1Map>::value_type threshold,
	     Compare compare)
{
    typedef typename vigragraph::graph_traits<Graph>::vertex_iterator src_scanner;
    typedef typename vigragraph::graph_traits<Graph2>::vertex_iterator dest_scanner;
    typedef typename vigragraph::graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename vigragraph::property_traits<T1Map>::value_type T1;
    typedef typename vigragraph::property_traits<T2Map>::value_type T2;

    src_scanner srcit, srcend;
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(Gsrc);

    dest_scanner destit, destend;
    vigragraph::tie(destit, destend) = vigragraph::vertices(Gdest);

    neighbor_iterator nbit, nbend;
#ifdef TEST_ONLY    
    size_t counter=0;
    size_t counter1=0;
    T1 dummysum=0;
#endif
    for (; srcit != srcend; ++srcit, ++destit) {
#ifdef TEST_ONLY    
      ++counter1;
#endif
	  const T1 refval = src[*srcit];

	  if (!compare(refval, threshold))
	      continue;
	  
	  vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*srcit, Gsrc);
	  
	  bool local_extremum = true;
	  int j=0;
	  for (;nbit != nbend; ++nbit) {
	    ++j;
#ifdef TEST_ONLY    
	    ++counter;
#endif
	    if (!compare(refval, src[*nbit])) {
	      local_extremum = false;
	      break;
	    }
	  }
	  if (local_extremum)
#ifdef TEST_ONLY    
	      dummysum += marker;
#else	 
  	      dest[*destit] = marker;
#endif

	  // TODO: case of no neighbors => local min/max per definition?
    }
#ifdef TEST_ONLY    
    std::cout << " counter=" << counter << " counter1=" << counter1 << std::endl;
    std::cout << " dummysum=" << ((int)dummysum) << std::endl;
#endif

}



  // Two-Iterator-Variant of the above: (vigra version of the BGL-like interface: vigra::ajdacent_vertices takes iterator, not *iterator argument.)
    template <class Graph, class Graph2, class T1Map, class T2Map, class Compare>
void
localMinMaxGraph2vigra(
	     Graph &Gsrc, 
 	     T1Map &src,
	     Graph2 &Gdest, 
 	     T2Map &dest,
	     typename vigragraph::property_traits<T2Map>::value_type marker, // unsigned int neighborhood,
	     typename vigragraph::property_traits<T1Map>::value_type threshold,
	     Compare compare)
{
    typedef typename vigragraph::graph_traits<Graph>::vertex_iterator src_scanner;
    typedef typename vigragraph::graph_traits<Graph2>::vertex_iterator dest_scanner;
    typedef typename vigragraph::graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename vigragraph::property_traits<T1Map>::value_type T1;
    typedef typename vigragraph::property_traits<T2Map>::value_type T2;

    src_scanner srcit, srcend;
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(Gsrc);

    dest_scanner destit, destend;
    vigragraph::tie(destit, destend) = vigragraph::vertices(Gdest);

    neighbor_iterator nbit, nbend;
#ifdef TEST_ONLY    
    size_t counter=0;
    size_t counter1=0;
    T1 dummysum=0;
#endif
    for (; srcit != srcend; ++srcit, ++destit) {
#ifdef TEST_ONLY    
      ++counter1;
#endif
      const T1 refval =  vigragraph::get(src, *srcit); // src[*srcit];

	  if (!compare(refval, threshold))
	      continue;
	  
#if 1
	  //vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*srcit, Gsrc);
	  vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices_at_iterator(srcit, Gsrc);

	  // std::pair<neighbor_iterator, neighbor_iterator> nbs = vigra::adjacent_vertices(srcit, Gsrc);
	  //	  neighbor_iterator &nbit(nbs.first), &nbend(nbs.second);
#else
	  // direct call:
	  nbit = Gsrc.get_neighbor_vertex_iterator(srcit);
	  nbend = Gsrc.get_neighbor_vertex_end_iterator(srcit);
#endif    
	  bool local_extremum = true;
	  int j=0;
	  for (;nbit != nbend; ++nbit) {
	    ++j;
#ifdef TEST_ONLY    
	    ++counter;
#endif
	    //	    if (!compare(refval, src[*nbit])) {
	    if (!compare(refval, vigragraph::get(src, *nbit))) {
	      local_extremum = false;
	      break;
	    }
	  }
	  if (local_extremum)
#ifdef TEST_ONLY    
	      dummysum += marker;
#else	 
	  vigragraph::put(dest, *destit, marker); //  dest[*destit] = marker;
#endif

	  // TODO: case of no neighbors => local min/max per definition?
    }
#ifdef TEST_ONLY    
    std::cout << " counter=" << counter << " counter1=" << counter1 << std::endl;
    std::cout << " dummysum=" << ((int)dummysum) << std::endl;
#endif

}


  // Two-Iterator-Variant of the above: (BGL-like interface)
    template <class Graph, class Graph2, class T1Map, class T2Map, class Compare>
void
localMinMaxGraph2boost(
	     Graph &Gsrc, 
 	     T1Map &src,
	     Graph2 &Gdest, 
 	     T2Map &dest,
	     typename vigragraph::property_traits<T2Map>::value_type marker, // unsigned int neighborhood,
	     typename vigragraph::property_traits<T1Map>::value_type threshold,
	     Compare compare)
{
    typedef typename vigragraph::graph_traits<Graph>::vertex_iterator src_scanner;
    typedef typename vigragraph::graph_traits<Graph2>::vertex_iterator dest_scanner;
    typedef typename vigragraph::graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename vigragraph::property_traits<T1Map>::value_type T1;
    typedef typename vigragraph::property_traits<T2Map>::value_type T2;

    src_scanner srcit, srcend;
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(Gsrc);

    dest_scanner destit, destend;
    vigragraph::tie(destit, destend) = vigragraph::vertices(Gdest);

    neighbor_iterator nbit, nbend;
#ifdef TEST_ONLY    
    size_t counter=0;
    size_t counter1=0;
    T1 dummysum=0;
#endif
    for (; srcit != srcend; ++srcit, ++destit) {
#ifdef TEST_ONLY    
      ++counter1;
#endif
      const T1 refval =  vigragraph::get(src, *srcit); // src[*srcit];

	  if (!compare(refval, threshold))
	      continue;
	  
	  vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*srcit, Gsrc);
	  bool local_extremum = true;
	  int j=0;
	  for (;nbit != nbend; ++nbit) {
	    ++j;
#ifdef TEST_ONLY    
	    ++counter;
#endif
	    //	    if (!compare(refval, src[*nbit])) {
	    if (!compare(refval, vigragraph::get(src, *nbit))) {
	      local_extremum = false;
	      break;
	    }
	  }
	  if (local_extremum)
#ifdef TEST_ONLY    
	      dummysum += marker;
#else	 
	  vigragraph::put(dest, *destit, marker); //  dest[*destit] = marker;
#endif

	  // TODO: case of no neighbors => local min/max per definition?
    }
#ifdef TEST_ONLY    
    std::cout << " counter=" << counter << " counter1=" << counter1 << std::endl;
    std::cout << " dummysum=" << ((int)dummysum) << std::endl;
#endif

}





  // Attempt without LValue propmaps, using only the free functions
  // to access ReadablePropertyMap (input) and WritablePropertyMap (label)
    template <class Graph, class T1Map, class T2Map, class Compare>
void
localMinMaxGraph3boost(
	     Graph const &G, 
 	     T1Map const &src,
 	     T2Map &dest,
	     typename vigragraph::property_traits<T2Map>::value_type marker, // unsigned int neighborhood,
	     typename vigragraph::property_traits<T1Map>::value_type threshold,
	     Compare const &compare)
{
    typedef typename vigragraph::graph_traits<Graph>::vertex_iterator graph_scanner;
    typedef typename vigragraph::graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename vigragraph::property_traits<T1Map>::value_type T1;
    typedef typename vigragraph::property_traits<T2Map>::value_type T2;

    graph_scanner srcit, srcend;
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(G);

    neighbor_iterator nbit, nbend;
#ifdef TEST_ONLY    
    size_t counter=0;
    size_t counter1=0;
    T1 dummysum=0;
#endif
    for (; srcit != srcend; ++srcit) {
#ifdef TEST_ONLY    
      ++counter1;
#endif
      const T1 refval = vigragraph::get(src, *srcit);

      if (!compare(refval, threshold))
	  continue;
	  
      // MAIN PROBLEM WITH BOOST INTERFACE:
      // adjacent_vertices is called with a vertex_descriptor,
      // not the iterator which possibly has more state!
      // -> potentially expensive to reconstruct iterator!
      vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*srcit, G);
	  
      bool local_extremum = true;
      int j=0;
      for (;nbit != nbend; ++nbit) {
	  ++j;
#ifdef TEST_ONLY    
	  ++counter;
#endif
	  if (!compare(refval, vigragraph::get(src, *nbit))) {
	      local_extremum = false;
	      break;
	  }
      }
	  if (local_extremum)
#ifdef TEST_ONLY    
	      dummysum += marker;
#else	 
	  vigragraph::put(dest, *srcit, marker);
#endif

	  // TODO: case of no neighbors => local min/max per definition?
    }
#ifdef TEST_ONLY    
    std::cout << " counter=" << counter << " counter1=" << counter1 << std::endl;
    std::cout << " dummysum=" << ((int)dummysum) << std::endl;
#endif

}




// not BGL-compatible: calls vigra::adjacent_vertices with vertex_iterator directly
// to avoid unnecessary conversion vertex_iterator (SSOI/CI) -> vertex_descriptor (coords) -> vertex_iterator (begin())+coord-offset
// TODO: also try to avoid "tie(...,...)" and use an iterator with hasMore() interface!
    template <class Graph, class T1Map, class T2Map, class Compare>
void
localMinMaxGraph3vigra(
	     Graph &G, 
 	     T1Map &src,
 	     T2Map &dest,
	     typename vigragraph::property_traits<T2Map>::value_type marker, // unsigned int neighborhood,
	     typename vigragraph::property_traits<T1Map>::value_type threshold,
	     Compare compare)
{
    typedef typename vigragraph::graph_traits<Graph>::vertex_iterator graph_scanner;
    typedef typename vigragraph::graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename vigragraph::property_traits<T1Map>::value_type T1;
    typedef typename vigragraph::property_traits<T2Map>::value_type T2;

    graph_scanner srcit, srcend;
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(G);

    neighbor_iterator nbit, nbend;
#ifdef TEST_ONLY    
    size_t counter=0;
    size_t counter1=0;
    T1 dummysum=0;
#endif
    for (; srcit != srcend; ++srcit) {
#ifdef TEST_ONLY    
      ++counter1;
#endif
      const T1 refval = vigragraph::get(src, *srcit);

      if (!compare(refval, threshold))
	  continue;
	  
      // MAIN PROBLEM WITH BOOST INTERFACE:
      // adjacent_vertices is called with a vertex_descriptor,
      // not the iterator which possibly has more state!
      // -> potentially expensive to reconstruct iterator!
      vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices_at_iterator(srcit, G);
	  
      bool local_extremum = true;
      int j=0;
      for (;nbit != nbend; ++nbit) {
	  ++j;
#ifdef TEST_ONLY    
	  ++counter;
#endif
	  if (!compare(refval, vigragraph::get(src, *nbit))) {
	      local_extremum = false;
	      break;
	  }
      }
	  if (local_extremum)
#ifdef TEST_ONLY    
	      dummysum += marker;
#else	
	      vigragraph::put(dest, *srcit, marker);
#endif

	  // TODO: case of no neighbors => local min/max per definition?
    }
#ifdef TEST_ONLY    
    std::cout << " counter=" << counter << " counter1=" << counter1 << std::endl;
    std::cout << " dummysum=" << ((int)dummysum) << std::endl;
#endif

}








template <class Graph, class T1Map, class T2Map, class TmpLabelMap, class Compare, class Equality>
void
extendedLocalMinMaxGraph(
	     Graph const &G, 
 	     T1Map const &src,
 	     T2Map &dest,
	     TmpLabelMap &tmpLabels,
	     typename vigragraph::property_traits<T2Map>::value_type marker, 
	     typename vigragraph::property_traits<T1Map>::value_type threshold,
	     Compare const &compare,
	     Equality const &equal)
{
    typedef typename vigragraph::graph_traits<Graph>::vertex_iterator graph_scanner;
    typedef typename vigragraph::graph_traits<Graph>::adjacency_iterator neighbor_iterator;

    typedef typename vigragraph::property_traits<T1Map>::value_type T1;
    typedef typename vigragraph::property_traits<T2Map>::value_type T2;


    int number_of_regions = labelGraph(G, src, tmpLabels, equal);
    // std::cout << "NUMBER OF REGIONS: " << number_of_regions << std::endl;

    // assume that a region is a extremum until the opposite is proved
    std::vector<unsigned char> isExtremum(number_of_regions + 1, (unsigned char)1);

    graph_scanner srcit, srcend;
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(G);
    neighbor_iterator nbit, nbend;
    for (; srcit != srcend; ++srcit) {
      const T1 refval = vigragraph::get(src, *srcit);

      const int lab = vigragraph::get(tmpLabels, *srcit);

      if (isExtremum[lab] == 0)
	  continue;

      if (!compare(refval, threshold)) {
	  // mark all regions that don't exceed the threshold as non-extremum
	  isExtremum[lab] = 0;
	  continue;	  
      }

      vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*srcit, G);
	  
      for (;nbit != nbend; ++nbit) {
	  const int nblab = vigragraph::get(tmpLabels, *nbit);
	  //	  if ((lab != nblab) && (!compare(refval, vigragraph::get(src, *nbit)))) {
	  if ((lab != nblab) && (compare(vigragraph::get(src, *nbit), refval))) {
	      isExtremum[lab] = 0;
	      break;
	  }
      }
    }
    
    // second pass: read out and mark extrema
    vigragraph::tie(srcit, srcend) = vigragraph::vertices(G);
    for (; srcit != srcend; ++srcit) {
	if (isExtremum[vigragraph::get(tmpLabels, *srcit)])
	    vigragraph::put(dest, *srcit, marker);
    }
}




} // namespace vigra



#endif // VIGRA_MULTI_LOCALMINMAX_HXX
