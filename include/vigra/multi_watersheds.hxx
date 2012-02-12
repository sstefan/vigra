#ifndef VIGRA_MULTI_WATERSHEDS_HXX
#define VIGRA_MULTI_WATERSHEDS_HXX

#include <vigra/voxelneighborhood.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/multi_localminmax.hxx>
#include <vigra/multi_labelgraph.hxx>
#include <vigra/multi_seededregiongrowing.hxx>
#include <vigra/watersheds.hxx>
#include <vigra/graphs.hxx>

namespace vigra
{


    template <class GRAPH, 
	      class SrcMap, 
	      class DestMap,
	      class TmpMap
	      >
    unsigned int
    generateWatershedSeeds_graph(GRAPH const &graph,
				 SrcMap const &src,
				 DestMap &dest,
				 TmpMap &tmp,
				 SeedOptions const & options = SeedOptions())
    {
	using namespace functor;
    
	typedef typename SrcMap::value_type value_type;
	typedef typename SrcMap::value_type marker_type;

	typedef typename GRAPH::vertex_descriptor VD;
	typedef typename GRAPH::vertex_iterator VI;
	typedef typename GRAPH::adjacency_iterator AI;


	vigra_precondition(options.mini != SeedOptions::LevelSets || 
			   options.thresholdIsValid<value_type>(),
			   "generateWatershedSeeds_graph(): SeedOptions.levelSets() must be specified with threshold.");
    
	if(options.mini == SeedOptions::LevelSets)
	    {
		VI ir, irend;
		vigragraph::tie(ir, irend) = vigragraph::vertices(graph);
		for (;ir!=irend; ++ir)
		    vigragraph::put(tmp, *ir, (vigragraph::get(src, *ir) <= options.thresh) ? 1 : 0);
	    }
	else
	    {
		// initialize temporary label map
		{
		    VI ir, irend;
		    vigragraph::tie(ir, irend) = vigragraph::vertices(graph);
		    for (;ir!=irend; ++ir)
			vigragraph::put(tmp, *ir, 0);
		}
		if (options.mini == SeedOptions::ExtendedMinima) {
		    extendedLocalMinMaxGraph(graph,
					     src, // =weights,
					     tmp, //=markers output,
					     dest, // use dest as temporary
					     marker_type(1), // marker
					     value_type(options.thresh), // threshold
					     std::less<value_type>(),
					     std::equal_to<value_type>());
		} else {
		    localMinMaxGraph3boost(graph,
					   src, // =weights,
					   tmp, //=markers output,
					   marker_type(1), // marker
					   value_type(options.thresh), // threshold
					   std::less<value_type>());
		}
	    
		if (0) { // TEST ONLY: write out seed labeling
		    VI i1, i1e;
		    vigragraph::tie(i1, i1e) = vigragraph::vertices(graph);
		    //char l='a';
		    int l=0;
		    for (;i1!=i1e; ++i1) {
			std::cout << " " << l;
			++l;
		    }
		    std::cout << std::endl;
		    vigragraph::tie(i1, i1e) = vigragraph::vertices(graph);
		    for (;i1!=i1e; ++i1) {
			std::cout << " " << vigragraph::get(tmp, *i1);
		    }
		    std::cout << std::endl;
		    vigragraph::tie(i1, i1e) = vigragraph::vertices(graph);
		    for (;i1!=i1e; ++i1) {
			std::cout << " " << vigragraph::get(src, *i1);
		    }
		    std::cout << std::endl;
		    std::cout << " thresh: " << options.thresh << std::endl;
		}
	    }
    
	const int numSeedRegions = 
	    labelGraphWithBackground(graph,
				     tmp,
				     dest,
				     0);
	if (0) {
	    std::cout << "Seeds after cc-labeling:" << std::endl;
	    VI i1, i1e;
	    vigragraph::tie(i1, i1e) = vigragraph::vertices(graph);
	    char l='a';
	    for (;i1!=i1e; ++i1) {
		std::cout << " " << l;
		++l;
	    }
	    std::cout << std::endl;
	    vigragraph::tie(i1, i1e) = vigragraph::vertices(graph);
	    for (;i1!=i1e; ++i1) {
		std::cout << " " << vigragraph::get(dest, *i1);
	    }
	    std::cout << std::endl;
	}
	// std::cout << " FOUND " << numSeedRegions << " SEEDS for watershed." << std::endl;
	return numSeedRegions;
    }


    template <class GRAPH, 
	      class SrcMap, 
	      class DestMap,
	      class TmpRegionsMap
	      >
    unsigned int
    watershedsRegionGrowing_graph(GRAPH const &graph,
				  SrcMap const &src,
				  DestMap &dest,
				  TmpRegionsMap &regions,
				  WatershedOptions const & options = WatershedOptions())
    {
	typedef typename SrcMap::value_type ValueType; 
	typedef typename DestMap::value_type LabelType; 
    
	unsigned int max_region_label = 0;
    
	if(options.seed_options.mini != SeedOptions::Unspecified)
	    {
		// we are supposed to compute seeds
		max_region_label = 
		    generateWatershedSeeds_graph(graph,
						 src, 
						 dest, 
						 regions, // used as temporary labeling array
						 options.seed_options);
	    }
    
	if(options.biased_label != 0)
	    {
		// create a statistics functor for biased region growing
		detail::BiasedWatershedStatistics<ValueType, LabelType> 
		    regionstats(options.biased_label, options.bias);

		// perform region growing, starting from the seeds computed above
		if(options.bucket_count == 0)
		    {
			max_region_label = 
			    seededRegionGrowing_graph(graph,
						      src,
						      dest,
						      dest,
						      regions,
						      regionstats, 
						      options.terminate,  // KeepContours?
						      options.max_cost);
		    }
		else
		    {
			max_region_label = 
			    fastSeededRegionGrowing_graph(graph,
							  src,  // source weights
							  dest, // seeds
							  dest, // destination
							  regionstats, 
							  options.terminate, 
							  options.max_cost, 
							  options.bucket_count);
		    }
	    }
	else
	    {
		// create a statistics functor for region growing
		detail::WatershedStatistics<ValueType, LabelType> regionstats;

		// perform region growing, starting from the seeds computed above
		if(options.bucket_count == 0)
		    {
			max_region_label =
			    seededRegionGrowing_graph(graph,
						      src,
						      dest,
						      dest,
						      regions,
						      regionstats, 
						      options.terminate,
						      options.max_cost);
		    }
		else
		    {
			max_region_label = 
			    fastSeededRegionGrowing_graph(graph,
							  src,  // source weights
							  dest, // seeds
							  dest, // destination
							  regionstats, 
							  options.terminate, 
							  options.max_cost, 
							  options.bucket_count);
		    }
	    }
    
	return max_region_label;
    }


}//namespace vigra

#endif /* VIGRA_MULTI_WATERSHEDS_HXX */
