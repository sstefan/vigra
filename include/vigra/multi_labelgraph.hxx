/************************************************************************/
/*                                                                      */
/*     Copyright 2006-2011 by Stefan Schmidt, F. Heinrich, B. Seppke,   */
/*                                                    Ullrich Koethe    */
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

#ifndef VIGRA_MULTI_LABELGRAPH_HXX
#define VIGRA_MULTI_LABELGRAPH_HXX


#include "graphs.hxx"
#include "voxelneighborhood.hxx"
#include "multi_array.hxx"
#include "union_find.hxx"

namespace vigra{




template <class GRAPH, 
	  class SrcMap, 
	  class DestMap,
	  class EqualityFunctor>
unsigned int labelGraph(GRAPH const &graph,
			SrcMap const &src,
			DestMap &dest,
			EqualityFunctor equal)
{
    typedef typename SrcMap::value_type ValueType;
    typedef typename DestMap::value_type LabelType;

    typedef typename GRAPH::vertex_descriptor VD;
    typedef typename GRAPH::vertex_iterator VI;
    typedef typename GRAPH::adjacency_iterator AI;

    // Use the ID map of the graph here to access SOI.
    typedef typename vigragraph::property_map<GRAPH, vigragraph::vertex_index_t>::type NodeIDMap;
    NodeIDMap nodeIDMap = vigragraph::get(vigragraph::vertex_index, graph);
                        
    // temporary image to store region labels
    detail::UnionFindArray<LabelType>  label;

    // pass 1: scan graph to find connected components

    // Each component will be represented by a tree of nodes. Each
    // node contains the scan order address of its parent in the
    // tree.  

    // In order for pass 2 to work correctly, the parent must
    // always have a smaller scan order address than the child.
    // Therefore, we can merge trees only at their roots, because the
    // root of the combined tree must have the smallest scan order
    // address among all the tree's pixels/ nodes.  The root of each
    // tree is distinguished by pointing to itself (it contains its
    // own scan order address). This condition is enforced whenever a
    // new region is found or two regions are merged

    // iterate source and destination graph:
    VI vi, viend;
    vigragraph::tie(vi, viend) = vigragraph::vertices(graph);
    for (;vi!=viend; ++vi) {
	ValueType curval = vigragraph::get(src, *vi);

#if 0  // the only difference to the code below... TODO: unify! (same in original labelvolume.hxx)
	if(equal(curval, backgroundValue))
	    {
		vigragraph::put(dest, *vi, label[0]);
		continue;
	    }
#endif

	LabelType currentLabel = label.nextFreeLabel();

	AI nbit, nbend;
	vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*vi, graph);
	for (;nbit != nbend; ++nbit) {
	    // check if "causal" neighbor by verifying if vertex ID is larger
	    // FIXME: The way this check is implemented may well be very expensive!
	    if (vigragraph::get(nodeIDMap, *nbit) < vigragraph::get(nodeIDMap, *vi)) {
		if(equal(curval, vigragraph::get(src, *nbit)))
		    {
			currentLabel = label.makeUnion(label[vigragraph::get(dest,*nbit)], currentLabel);
		    }
	    }
	}

	vigragraph::put(dest, *vi, label.finalizeLabel(currentLabel));
    }
    
    LabelType count = label.makeContiguous();

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    
    vigragraph::tie(vi, viend) = vigragraph::vertices(graph);
    for (;vi!=viend; ++vi) 
	vigragraph::put(dest, *vi, label[vigragraph::get(dest, *vi)]);
    return count;
}



template <class GRAPH, 
	  class SrcMap, 
	  class DestMap>
unsigned int labelGraph(GRAPH const &graph,
			SrcMap const &src,
			DestMap  &dest)
{
    return labelGraph(graph,
		      src,
		      dest,
		      std::equal_to<typename SrcMap::value_type>());
}

/********************************************************/
/*                                                      */
/*               labelGraphWithBackground              */
/*                                                      */
/********************************************************/

/** \brief Find the connected components of a segmented graph,
     excluding the background from labeling.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
    ... TODO ...
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
    ... TODO ...
    }
    \endcode

    Connected components are defined as regions with uniform node
    values. Thus, the input graph property map either must be
    equality comparable (first form), or an EqualityFunctor must be
    provided that realizes the desired predicate (second form). All
    nodex equal to the given '<TT>background_value</TT>' are ignored
    when determining connected components and are set to label 0 in the
    destination property map.

    The destination's value type should be large enough to hold the
    labels without overflow. Region numbers will be a consecutive
    sequence starting with one and ending with the region number
    returned by the function (inclusive).

    Return:  the number of regions found (= largest region label)

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_labelgraph.hxx\><br>
    Namespace: vigra

    \code
    ... TODO ...
    \endcode

*/

template <class GRAPH, 
	  class SrcMap, 
	  class DestMap,
	  class ValueType,
	  class EqualityFunctor>
unsigned int labelGraphWithBackground(GRAPH const &graph,
				       SrcMap const &src,
				       DestMap &dest,
                                       ValueType backgroundValue,
				       EqualityFunctor equal)
{
    typedef typename DestMap::value_type LabelType;

    typedef typename GRAPH::vertex_descriptor VD;
    typedef typename GRAPH::vertex_iterator VI;
    typedef typename GRAPH::adjacency_iterator AI;

    // Use the ID map of the graph here to access SOI.
    typedef typename vigragraph::property_map<GRAPH, vigragraph::vertex_index_t>::type NodeIDMap;
    NodeIDMap nodeIDMap = vigragraph::get(vigragraph::vertex_index, graph);
                        
    // temporary image to store region labels
    detail::UnionFindArray<LabelType>  label;

    // pass 1: scan graph to find connected components

    // Each component will be represented by a tree of nodes. Each
    // node contains the scan order address of its parent in the
    // tree.  

    // In order for pass 2 to work correctly, the parent must
    // always have a smaller scan order address than the child.
    // Therefore, we can merge trees only at their roots, because the
    // root of the combined tree must have the smallest scan order
    // address among all the tree's pixels/ nodes.  The root of each
    // tree is distinguished by pointing to itself (it contains its
    // own scan order address). This condition is enforced whenever a
    // new region is found or two regions are merged

    // iterate source and destination graph:
    VI vi, viend;
    vigragraph::tie(vi, viend) = vigragraph::vertices(graph);
    for (;vi!=viend; ++vi) {
	ValueType curval = vigragraph::get(src, *vi);
	if(equal(curval, backgroundValue))
	    {
		vigragraph::put(dest, *vi, label[0]);
		continue;
	    }

	LabelType currentLabel = label.nextFreeLabel();

	AI nbit, nbend;
	vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*vi, graph);
	for (;nbit != nbend; ++nbit) {
	    // check if "causal" neighbor by verifying if vertex ID is larger
	    // FIXME: The way this check is implemented may well be very expensive!
	    if (vigragraph::get(nodeIDMap, *nbit) < vigragraph::get(nodeIDMap, *vi)) {
		if(equal(curval, vigragraph::get(src, *nbit)))
		    {
			currentLabel = label.makeUnion(label[vigragraph::get(dest,*nbit)], currentLabel);
		    }
	    }
	}

	vigragraph::put(dest, *vi, label.finalizeLabel(currentLabel));
    }
    
    LabelType count = label.makeContiguous();

    // pass 2: assign one label to each region (tree)
    // so that labels form a consecutive sequence 1, 2, ...
    
    vigragraph::tie(vi, viend) = vigragraph::vertices(graph);
    for (;vi!=viend; ++vi) 
	vigragraph::put(dest, *vi, label[vigragraph::get(dest, *vi)]);
    return count;
}



template <class GRAPH, 
	  class SrcMap, 
	  class DestMap,
	  class ValueType>
unsigned int labelGraphWithBackground(GRAPH const &graph,
				       SrcMap const &src,
				       DestMap &dest,
                                       ValueType backgroundValue)
{
    return labelGraphWithBackground(graph,
				     src,
				     dest,
                                     backgroundValue,
				     std::equal_to<typename SrcMap::value_type>());
}



//@}

} //end of namespace vigra

#endif //VIGRA_MULTI_LABELGRAPH_HXX
