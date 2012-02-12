/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2012 by Kasim Terzic, Christian-Dennis Rahn,      */
/*                        Stefan Schmidt, and Ullrich Koethe            */
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

#ifndef VIGRA_MULTI_SEEDEDREGIONGROWING_HXX
#define VIGRA_MULTI_SEEDEDREGIONGROWING_HXX

#include <vector>
#include <stack>
#include <queue>
#include <vigra/utilities.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <vigra/voxelneighborhood.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/metaprogramming.hxx>
#include <vigra/graphs.hxx>


namespace vigra {

namespace detail {

    // FIXME: in general graphs, it doesn't make too much sense to compute the distance via indices.
    // We should rather use some special (BGL?) distance function for this.
    // For now, we overload "dot" for scalar types, such that the code below also works for non-grid-graphs.
    // TODO: The tie-breaker code should be refactored altogether. The vertex index / coordinates should be 
    // handed over to some cost functor and tie-braking may be handled there...
     template<class T>
     typename PromoteTraits<T, T>::Promote
     dot(T x1, typename If<typename NumericTraits<T>::isScalar, T, void>::type x2)
     {
 	return x1*x2;
     }



    template <class GRAPH,
	      class COST, 
	      class Location_type,
	      class Nearest_type,
	      class SrcType,
	      class LabelType>
    class SeedRgQueueObjectGridGraph
    {
    public:
	typedef SeedRgQueueObjectGridGraph<GRAPH,COST,Location_type,Nearest_type,SrcType,LabelType> queueObject_type;

	Location_type location_;
	Nearest_type nearest_;
	COST cost_;
	int count_;
	int label_;
	int dist_;

	static inline 
	int dist(Location_type current, Nearest_type nearest) {
	    return 	dot(current-nearest, current-nearest);
	}

	SeedRgQueueObjectGridGraph()
	//    	: location_(0), nearest_(0), cost_(0), count_(0), label_(0)
	{
#if 1
	    // this seems to be much faster than the initializer list above with -O2.... no difference with -O3
// 	    location_ = Location_type(0);
// 	    nearest_ = Nearest_type(0);
	    cost_ = 0;
	    count_ = 0;
	    label_ = 0;
#endif
	}

	SeedRgQueueObjectGridGraph(Location_type const & location, Nearest_type const & nearest,
				   COST const & cost, int const & count, int const & label)
	    : location_(location), nearest_(nearest),
	      cost_(cost), count_(count), label_(label)
	{
	    dist_ = dist(location_, nearest_);
	}

	void set(Location_type const & location, Nearest_type const & nearest,
		 COST const & cost, int const & count, int const & label)
	{
	    location_ = location;
	    nearest_ = nearest;
	    cost_ = cost;
	    count_ = count;
	    label_ = label;

	    dist_ = dist(location_, nearest_);
	}


	struct Compare
	{
	    // must implement > since priority_queue looks for largest element
	    bool operator()(SeedRgQueueObjectGridGraph const & l,
			    SeedRgQueueObjectGridGraph const & r) const
	    {
		if(r.cost_ == l.cost_)
		    {
			if(r.dist_ == l.dist_) return r.count_ < l.count_;

			return r.dist_ < l.dist_;
		    }

		return r.cost_ < l.cost_;
	    }
	    bool operator()(SeedRgQueueObjectGridGraph const * l,
			    SeedRgQueueObjectGridGraph const * r) const
	    {
		if(r->cost_ == l->cost_)
		    {
			if(r->dist_ == l->dist_) return r->count_ < l->count_;

			return r->dist_ < l->dist_;
		    }

		return r->cost_ < l->cost_;
	    }
	};




	struct Allocator
	{
	    ~Allocator()
	    {
		while(!freelist_.empty())
		    {
			delete freelist_.top();
			freelist_.pop();
		    }
	    }



	    SeedRgQueueObjectGridGraph * create(Location_type const & location, Nearest_type const & nearest,
						COST const & cost, int const & count, int const & label)
	    {
		if(!freelist_.empty())
		    {
			SeedRgQueueObjectGridGraph * res = freelist_.top();
			freelist_.pop();
			res->set(location, nearest, cost, count, label);
			return res;
		    }

		return new SeedRgQueueObjectGridGraph(location, nearest, cost, count, label);
	    }


	    void dismiss(SeedRgQueueObjectGridGraph * p)
	    {
		freelist_.push(p);
	    }

	    std::stack<queueObject_type *> freelist_;
	};



	class QueueObjectContext
	{
	public:
	    typedef queueObject_type QueueObject;
	    typedef Nearest_type Origin_type;
	    typedef typename GRAPH::vertex_iterator label_iterator;
	    typedef typename GRAPH::adjacency_iterator label_nb_iterator;
	    typedef typename GRAPH::vertex_iterator source_iterator;
	    typedef typename GRAPH::adjacency_iterator source_nb_iterator;

	    QueueObjectContext(GRAPH const &graph, SrcType const &src, LabelType &labels) :
		graph_(graph), src_(src), labels_(labels)
	    {}

	    inline int getLabel(const QueueObject &qo) const 
	    {
		return qo.label_;
	    }
	    inline COST getCost(const QueueObject &qo) const 
	    {
		return qo.cost_;
	    }
	    inline Origin_type getOrigin(const QueueObject &qo) const
	    {
		return qo.nearest_;
	    }

	    label_iterator getLabelImageIterator(const QueueObject &qo) const
	    {
		// TODO:
		// test a variant where data needed for this reconstruction (scanorderindex) is available in the QueueObject

		// reconstruct a full strided scan order iterator 
		// (at least any iterator providing the data pointer AND knowing the neighborhood at the point)
		// FIXME: probably very inefficient!
		return vigragraph::vertices(graph_).first + qo.location_;
		// 	    label_iterator tmp = vigragraph::vertices(graph_).first;
		// 	    tmp += qo.location_; // CHECK: not a valid expression for BGL iterators?
		// 	    return tmp;
	    }

	    inline
	    std::pair<label_nb_iterator,label_nb_iterator> 
	    label_neighbor_iteratorpair(const QueueObject &qo, label_iterator& label_iterator) const
	    {
		return vigragraph::adjacent_vertices(qo.location_, graph_);
	    }


	    source_iterator getSrcImageIterator(const QueueObject &qo) const
	    {
		return vigragraph::vertices(graph_).first + qo.location_;
	    }

	    inline 
	    std::pair<source_nb_iterator,source_nb_iterator>
	    source_neighbor_iteratorpair(const QueueObject &qo, source_iterator& source_iterator) const
	    {
		return vigragraph::adjacent_vertices(qo.location_, graph_);
	    }


	protected:
	    // e.g. store GridGraph reference etc. here
	    //   const 
	    GRAPH const &graph_;
	    SrcType const &src_;
	    LabelType &labels_;
	};
    };




} // namespace detail









    //! generic SRG algorithm for BGL style graphs
    // FIXME: Needs cleanup and more testing!
    template <class GRAPH, 
	      class SrcMap, 
	      class SeedMap, 
	      class DestMap, 
	      class TmpRegionsMap, 
	      class RegionStatisticsArray>
    typename DestMap::value_type
    seededRegionGrowing_graph(GRAPH const &graph,
			      SrcMap const &src,
			      SeedMap const &seed,
			      DestMap &dest,
			      TmpRegionsMap &regions,
			      RegionStatisticsArray &stats, 
			      SRGType srgType,
			      double max_cost = 0)
{
    // BOOST_CONCEPT_ASSERT((AdjacencyGraphConcept<GRAPH>));

    typedef typename GRAPH::vertex_descriptor VD;
    typedef typename GRAPH::vertex_iterator VI;
    typedef typename GRAPH::adjacency_iterator AI;


    int count = 0;
    typedef typename RegionStatisticsArray::value_type RegionStatistics;
    typedef typename PromoteTraits<typename RegionStatistics::cost_type, double>::Promote CostType;

    typedef detail::SeedRgQueueObjectGridGraph<GRAPH, CostType, VD, VD, SrcMap, TmpRegionsMap> Node;
    typedef typename Node::QueueObjectContext QueueObjectContext;

    typename Node::Allocator allocator;
    QueueObjectContext queueObjectContext(graph, src, regions);

    typedef std::priority_queue< Node *,
	std::vector<Node *>,
	typename Node::Compare >  SeedRgNodeHeap;

    
    //initImageBorder(destImageRange(regions), 1, SRGWatershedLabel);
    // @ask FIXME: border handling in this way really required?
    //  (commented out to be able to use stridedscanorderiterators)
    //    initMultiArrayBorder(destMultiArrayRange(regions), 1, SRGWatershedLabel); 
    
    // copy seed image to our local label image:
    VI ir, irend;
    vigragraph::tie(ir, irend) = vigragraph::vertices(graph);
    for (;ir!=irend; ++ir)
	vigragraph::put(regions, *ir, vigragraph::get(seed, *ir));

    // allocate and init memory for the results
    SeedRgNodeHeap pheap;
    int labelOfNeighbor;
    int maxRegionLabel = 0;

    // Joint iteration over source image and label image 
    // using strided scan order iterators.
    // (Could also use some joint iterator over both the input image and the label image here)
    VI is, isend;
    vigragraph::tie(is, isend) = vigragraph::vertices(graph);
    for (; is!=isend; ++is) 
	{
	    if(vigragraph::get(regions, *is) == 0) // unlabeled point
		{
		    AI nbit, nbend;
		    vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*is, graph);

		    // find candidate pixels for growing and fill heap
		    int i=0;
		    for (;nbit != nbend; ++nbit)
			{
			    ++i;
			    labelOfNeighbor = vigragraph::get(regions, *nbit);
			    if(labelOfNeighbor > 0) // is neighbor to a labeled point -> push candidate on queue
				{
				    CostType cost = stats[labelOfNeighbor].cost(vigragraph::get(src, *is));
				    
				    Node * node =
					allocator.create(*is,
							 *nbit,
							 cost, count++, labelOfNeighbor);
				    pheap.push(node);
				}
			    // NOTE: This might produce several queue entries for the same node, in case it is neighbor to
			    //   multiple labeled points (but eventually with different costs) 
			    //   ... this is handled in the pop routine, by discarding already labeled points there.
			}
		}
	    else {
		vigra_precondition((typename DestMap::value_type)vigragraph::get(regions, *is) <= stats.maxRegionLabel(),
				   "seededRegionGrowing_graph(): Largest label exceeds size of RegionStatisticsArray.");
                if(maxRegionLabel < vigragraph::get(regions, *is))
                    maxRegionLabel = vigragraph::get(regions, *is);
		
	    }
	}



    typename QueueObjectContext::label_nb_iterator 
	lnbit, lnbend;

    typename QueueObjectContext::source_nb_iterator 
	snbit, snbend;


    // perform region growing
    while(pheap.size() != 0)
	{
	    Node * node = pheap.top();
	    pheap.pop();

	    typename QueueObjectContext::Origin_type nearest = queueObjectContext.getOrigin(*node);
	    int lab = queueObjectContext.getLabel(*node);
	    CostType cost = queueObjectContext.getCost(*node);

	    //	std::cerr << " pop lab=" << lab << " cost= " << cost << std::endl;	

	    // because we free the node next, we need a copy currently.
	    // alternative: pre-allocate the iterators required below
	    Node tmpnode(*node);

	    allocator.dismiss(node);

	    // @ask: couldn't we have just eliminated inserting it instead?
	    if((srgType & StopAtThreshold) != 0 && cost > max_cost) {
		std::cerr << "Stopping at threshold " << cost << std::endl; // FIXME: remove
		break;
	    }


	    typename QueueObjectContext::label_iterator irx = queueObjectContext.getLabelImageIterator(tmpnode);
	    typename QueueObjectContext::source_iterator isx = queueObjectContext.getSrcImageIterator(tmpnode); // also required in source image

	    if(vigragraph::get(regions,*irx)) // already labelled region / watershed?
		continue;

	    if((srgType & KeepContours) != 0)
		{
		    // iterate the respective neighbors of this position
		    // in the label image.
		    vigragraph::tie(lnbit, lnbend) = queueObjectContext.label_neighbor_iteratorpair(tmpnode, irx);

		    for(; lnbit != lnbend; ++lnbit)
			{
			    labelOfNeighbor = vigragraph::get(regions, *lnbit);
			    if((labelOfNeighbor>0) && (labelOfNeighbor != lab))
				{
				    // one of the current's points neighbors has already been discovered from a different basin.
				    // hence the current point must be a watershed:
				    lab = SRGWatershedLabel;
				    break;
				}
			}
		}

	    vigragraph::put(regions, *irx, lab);

	    if((srgType & KeepContours) == 0 || lab > 0)
		{
		    // update statistics
		    stats[vigragraph::get(regions, *irx)](vigragraph::get(src, *isx)); 
	    
		    // search neighborhood
		    // second pass: find new candidate pixels

		    // iterate the respective neighbors of this position
		    // in the label image AND the source image
		    vigragraph::tie(lnbit, lnbend) = queueObjectContext.label_neighbor_iteratorpair(tmpnode, irx);		    
		    vigragraph::tie(snbit, snbend) = queueObjectContext.source_neighbor_iteratorpair(tmpnode, isx);		    
	    
		    for(; lnbit != lnbend; ++lnbit,++snbit)
			{
			    if(vigragraph::get(regions, *lnbit) == 0)
				{
				    CostType cost = stats[lab].cost(vigragraph::get(src, *snbit));

				    Node * new_node =
					allocator.create(*lnbit, nearest, cost, count++, lab);
				    pheap.push(new_node);
				}
			}
		}
	}
    
    // free temporary memory
    while(pheap.size() != 0)
	{
	    allocator.dismiss(pheap.top());
	    pheap.pop();
	}

    // write result
    vigragraph::tie(is, isend) = vigragraph::vertices(graph);
    detail::UnlabelWatersheds transformer;
    for (; is!=isend; ++is) {
	vigragraph::put(dest, *is, transformer(vigragraph::get(regions, *is)));
    }
    return (typename DestMap::value_type) maxRegionLabel;
}



    //! Turbo SRG algorithm for BGL style graphs
    //  using Bucket Queue (turbo watershed algorithm, 
    //  modeled more closely after existing seededregiongrowing.hxx code)
    template <class GRAPH, 
	      class SrcMap, 
	      class SeedMap, 
	      class DestMap, 
	      class RegionStatisticsArray>
    typename DestMap::value_type
    fastSeededRegionGrowing_graph(GRAPH const &graph,
				  SrcMap const &src,
				  SeedMap const &seed,
				  DestMap &dest,
				  RegionStatisticsArray &stats, 
				  SRGType srgType,
				  double max_cost = 0,
				  std::ptrdiff_t bucket_count = 256)
{
    vigra_precondition((srgType & KeepContours) == 0,
       "fastSeededRegionGrowing(): the turbo algorithm doesn't support 'KeepContours', sorry.");

    typedef typename GRAPH::vertex_descriptor VD;
    typedef typename GRAPH::vertex_iterator VI;
    typedef typename GRAPH::adjacency_iterator AI;

    typedef typename DestMap::value_type LabelType;
    
    typedef typename RegionStatisticsArray::value_type RegionStatistics;
    typedef typename PromoteTraits<typename RegionStatistics::cost_type, double>::Promote CostType;

    typedef VD QueueObject;
    QueueObject tmpQueueObject;

    // define one queue per possible input level:
    std::vector<std::queue<QueueObject> > queues(bucket_count);	
    
    // copy seed image to our local label image:
    VI ir, irend;
    vigragraph::tie(ir, irend) = vigragraph::vertices(graph);
    for (;ir!=irend; ++ir)
	vigragraph::put(dest, *ir, vigragraph::get(seed, *ir));

    unsigned int maxRegionLabel = 0;

    VI is, isend;
    vigragraph::tie(is, isend) = vigragraph::vertices(graph);
    for (; is!=isend; ++is) {

	if(vigragraph::get(dest, *is) != 0) { // unlabeled point		
	    const LabelType currentLabel = vigragraph::get(dest, *is);
	    vigra_precondition(currentLabel <= stats.maxRegionLabel(),
			       "fastSeededRegionGrowing_graph(): Largest label exceeds size of RegionStatisticsArray.");
	    if(maxRegionLabel < currentLabel)
		maxRegionLabel = currentLabel;

	    AI nbit, nbend;
	    vigragraph::tie(nbit, nbend) = vigragraph::adjacent_vertices(*is, graph);

	    // find candidate pixels for growing and fill heap
	    int i=0;
	    for (;nbit != nbend; ++nbit) {
		++i;
		const LabelType labelOfNeighbor = vigragraph::get(dest, *nbit);
		if(labelOfNeighbor == 0) // neighbor is unlabeled -> push candidate on queue
		    {
			const CostType cost = stats[currentLabel].cost(vigragraph::get(src, *is));
			tmpQueueObject = QueueObject(*is);
			queues[cost].push(tmpQueueObject);

			break; // HERE: only push once if it has unlabeled neighbors.
			// they all have the same cost, because it's dependent on the current vertex only,
			// and not the neighbor.
		    }
	    }

	}
    }


    AI lnbit, lnbend;

    CostType priority = 0;
    // perform region growing
    while (priority < bucket_count) {
	while(!queues[priority].empty())
	{
	    const VD currentVertex = queues[priority].front();
	    queues[priority].pop();

	    int lab = vigragraph::get(dest, currentVertex);

	    // @ask: couldn't we have just eliminated inserting it instead?
	    if((srgType & StopAtThreshold) != 0 && priority > max_cost) {
		break;
	    }

	    // iterate neighbors
	    vigragraph::tie(lnbit, lnbend) =  vigragraph::adjacent_vertices(currentVertex, graph);
	    for(; lnbit != lnbend; ++lnbit) {
		const LabelType labelOfNeighbor = vigragraph::get(dest, *lnbit);
		if(labelOfNeighbor == 0) {

		    vigragraph::put(dest, *lnbit, lab); // label that neighbor
		    // and push it to the queue
		    const CostType cost = std::max(priority,
						   (CostType)stats[lab].cost(vigragraph::get(src, *lnbit)));
		    queues[cost].push(*lnbit);
		}
	    }	    
	} // while !empty
	queues[priority] = std::queue<QueueObject>(); // free processed bucket
	++priority;
    } // while true
    
    return (typename DestMap::value_type) maxRegionLabel;
}







} // namespace vigra

#endif // VIGRA_MULTI_SEEDEDREGIONGROWING_HXX

