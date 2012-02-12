//#define SIM_PQ_WATERSHED
#define SIM_TURBO_WATERSHED
// #define USE_COMPILED_TEST_DATA

#include <typeinfo>
#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <queue>

#include <unittest.hxx>

#include <vigra/timing.hxx>

// for multi-dim image experiments:
#include <vigra/multi_array.hxx>

// just an empty function with void* argument (in different compilation unit), used to prevent optimizing away code
#include "empty.hxx"


USETICTOC

using namespace std;
using namespace vigra;

const bool verbose = true;
const size_t repetitions = 1;
const size_t best_of_repetitions = 3;

typedef vigra::MultiArray<3, unsigned char> IMG;


#ifdef USE_COMPILED_TEST_DATA
//#include "log.cxx"
//#include "log_test.cxx"
#include "log.hxx"




struct TestData {
    typedef MultiArrayShape<3>::type Shp3D;
    typedef IMG image_type;
    image_type in;

    const size_t * const logdata;
    const size_t numentries; 
    const size_t *  iter;
    const size_t * const iter_end;
    const size_t maxval; // maximum occuring in test dataset 

    size_t counter;

    TestData() 
	: logdata(LogData::data), 
	  numentries(LogData::count), 
	  iter(logdata),
	  iter_end(logdata+3*LogData::count),
	  maxval(1792),
	  counter(0)
    {
	Shp3D shape(LogData::shape[0], LogData::shape[1], LogData::shape[2]);
	in = image_type(shape);
	// image "in" data will not actually be used...
	// we will rather simulate the queue access pattern based on the log data
    }
    
    bool getNext(bool &type, size_t &priority, size_t &id) {
	if (iter >= iter_end)
	    return false;
	++counter;
	type = *(iter++); // true means push, false means pop operation
	priority = *(iter++);
	id = *(iter++);
	return true;
    }

    void reset() {
	iter = logdata;
	counter = 0;
    }
};

#else

struct TestData {
    typedef MultiArrayShape<3>::type Shp3D;
    typedef IMG image_type;
    image_type in;
    vector<size_t> log;

    const size_t * logdata;
    size_t numentries; 
    const size_t * iter;
    const size_t * iter_end;
    size_t maxval;

    size_t counter;

    TestData() 
	: counter(0)
    {
	Shp3D shape;

	const string filename = "log.bin";
	std::ifstream inp(filename.c_str(), ifstream::in | ifstream::binary);
	unsigned int tmp;
	inp.read((char*)&tmp, sizeof(tmp)); shape[0] = tmp;
	inp.read((char*)&tmp, sizeof(tmp)); shape[1] = tmp;
	inp.read((char*)&tmp, sizeof(tmp)); shape[2] = tmp;
	inp.read((char*)&tmp, sizeof(tmp)); numentries = tmp;
	// cout << "Reading in " << numentries << " entries:" << endl;

	log.resize(3*numentries);
	size_t *plog = &log[0];
	maxval = 0;
	for (size_t i=0; i < numentries; ++i) {
	    inp.read((char*)&tmp, sizeof(tmp)); *(plog++) = tmp;
	    inp.read((char*)&tmp, sizeof(tmp)); *(plog++) = tmp;
	    if (tmp > maxval) 
		maxval = tmp;
	    inp.read((char*)&tmp, sizeof(tmp)); *(plog++) = tmp;
	}

	inp.close();

	logdata = iter = &log[0];
	iter_end = logdata+3*numentries;
	
	in = image_type(shape);
	// image "in" data will not actually be used...
	// we will rather simulate the queue access pattern based on the log data
    }
    
    bool getNext(bool &type, size_t &priority, size_t &id) {
	if (iter >= iter_end)
	    return false;
	++counter;
	type = *(iter++); // true means push, false means pop operation
	priority = *(iter++);
	id = *(iter++);
	return true;
    }

    void reset() {
	iter = logdata;
	counter = 0;
    }
};
#endif


// Functors to apply to heap objects to test representation conversion overhead:

//! No-Operation-Functor, as baseline to measure queue overhead itself
struct ConversionNOP {
    typedef size_t vertex_descriptor;
    typedef size_t queue_vertex_descriptor;
    typedef size_t algo_vertex_descriptor;

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	doNothing((void*)&in);
    }

    static inline void execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	//	vertex_descriptor ret = vertex_descriptor(in);
	doNothing((void*)&in);
    }

    static inline size_t convert(const IMG &img, const size_t &id)
    {
	return id;
    }

    static const string name() { return "NOP"; }
};



//! This class simulates the case that only the ScanOrderIndex is put on the queues.
//   and a full ScanOrderIterator is used by the Algorithm
struct SOI_OnQueue : ConversionNOP {
    typedef IMG::view_type view_type;

    // algorithm uses the full scridedscanorderiterator:
    typedef view_type::iterator algo_vertex_descriptor; // the StridedScanOrderIterator
    typedef view_type::iterator vertex_descriptor;

    typedef vigra::MultiArrayIndex queue_vertex_descriptor;

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	// on Pop, need to reconstruct coords
	// or even full SSOI to iterate neighborhood:
	out = static_cast<view_type>(img).begin() + in;
    }


    static inline void execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	// on Push, just extract SOI
	out = in.index();
    }


    static inline vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	vertex_descriptor tmp = static_cast<view_type>(img).begin();
	tmp += id;
	return tmp;
    }

    static const string name() { return "SOI_OnQueue"; }
};





//! This class simulates the case that only the coordinates tuple is put on the queues,
//  and a full scanorderiterator is used in the algorithm:
struct CoordsOnQueue : ConversionNOP {
    typedef IMG::view_type view_type;

    // uses the coordinates on the queue:
    typedef IMG::difference_type  queue_vertex_descriptor;
    // algorithm uses the full scridedscanorderiterator:
    typedef view_type::iterator algo_vertex_descriptor; // the StridedScanOrderIterator

    // the object processed by the algorithm, to be prepared (reconstructed from SOI) beforehand:
    typedef view_type::iterator vertex_descriptor; // the StridedScanOrderIterator

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    { 
	// convert coords to full scan order iterator:
	out = static_cast<view_type>(img).begin() + img.coordinateToScanOrderIndex(in);
    }

    static inline void execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    { 	// push executes for every unlabeled neighbor of a point.
	// we just store the (already available) coordinates here.
	out = in.point();
    }

    static inline vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	vertex_descriptor tmp = static_cast<view_type>(img).begin();
	tmp += id;
	return tmp;
    }

    static const string name() { return "CoordsOnQueue"; }
};









//! This class simulates the case that a full StridedScanOrderIterator
//  is put on the queue.
//   In addition to the coordinates tuple, it already has the 
//   neighborhood information and can iterate its neighbors,
//   and furthermore the linear index / ptr is also available directly.
struct IterOnQueue : ConversionNOP {
    typedef IMG::view_type view_type;
    typedef view_type::iterator vertex_descriptor; // the StridedScanOrderIterator
    typedef view_type::iterator queue_vertex_descriptor;
    typedef view_type::iterator algo_vertex_descriptor; 

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	// no special actions to take on pop. 
	// doNothing((void*)&in);
	//  however, we need to copy the object off the heap/queue:
	out = in;
    }

    static inline void execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	out = in; // FIXME: Can we avoid this copy/assignment as well?
	// FIXME: Instead of using the output argument, may directly pass reference to
	// input as return value?	
    }

    static inline vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	vertex_descriptor res = static_cast<view_type>(img).begin();
	res += id; 
	return res;
    }

    static const string name() { return "IterOnQueue"; }
};


//! This class simulates the case that a full StridedScanOrderIterator
//  is put on the queue.
//   In addition to the coordinates tuple, it already has the 
//   neighborhood information and can iterate its neighbors,
//   and furthermore the linear index / ptr is also available directly.
struct IterOnQueue_NoCopies : IterOnQueue {
    typedef IMG::view_type view_type;
    typedef view_type::iterator vertex_descriptor; // the StridedScanOrderIterator
    typedef view_type::iterator queue_vertex_descriptor;
    typedef view_type::iterator algo_vertex_descriptor; 

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	// no special actions to take on pop. 
	// doNothing((void*)&in);
	//  however, we need to copy the object off the heap/queue:
	//	out = in;
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
    }

    static inline void execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
    }

    static const string name() { return "IterOnQueue_NoCopies"; }
};




//! This class simulates the case that a full StridedScanOrderIterator
//  is put on the queue.
//   In addition to the coordinates tuple, it already has the 
//   neighborhood information and can iterate its neighbors,
//   and furthermore the linear index / ptr is also available directly.
struct IterOnQueue_ManyNOPs : IterOnQueue {
    typedef IMG::view_type view_type;
    typedef view_type::iterator vertex_descriptor; // the StridedScanOrderIterator
    typedef view_type::iterator queue_vertex_descriptor;
    typedef view_type::iterator algo_vertex_descriptor; 

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	// no special actions to take on pop. 
	// doNothing((void*)&in);
	//  however, we need to copy the object off the heap/queue:
	//	out = in;
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
    }

    static inline void execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
    }

    static const string name() { return "IterOnQueue_ManyNOPs"; }
};







template<class IDTYPE, int PAYLOADSIZE>
struct TurboQueueObject 
{
    IDTYPE id_;
    unsigned char data_[PAYLOADSIZE];
    TurboQueueObject(const IDTYPE &id)
	: id_(id) {
	for (size_t i=0; i<PAYLOADSIZE; i++)
	    data_[i] = 0;
    }

    // provide copy constructor and assignment operator
    // to more realistically simulate a heterogeneous object 
    // where simple memset/memcopy is not possible
    TurboQueueObject(const TurboQueueObject &copy_from) 
	: id_(copy_from.id_)
    {
	// FIXME: use template-code based loop unrolling here !?
	// or does the compiler optimize this loop already anyways?
	for (size_t i=0; i<PAYLOADSIZE; i++)
	    data_[i] = copy_from.data_[i];
    }

    inline TurboQueueObject& operator= (const TurboQueueObject &assign_from)
    {
	id_ = assign_from.id_;
	for (size_t i=0; i<PAYLOADSIZE; i++)
	    data_[i] = assign_from.data_[i];
	return *this;
    }
};

template<class WEIGHTTYPE, class IDTYPE, int PAYLOADSIZE>
struct HeapObject : TurboQueueObject<IDTYPE, PAYLOADSIZE>
{
    typedef TurboQueueObject<IDTYPE, PAYLOADSIZE> super_type;
    WEIGHTTYPE cost_;

    HeapObject(const WEIGHTTYPE &cost, const IDTYPE &id)
	: super_type(id), cost_(cost)  {}

    struct Compare
    {
        // must implement > since priority_queue looks for largest element
        bool operator()(HeapObject const & l,
                        HeapObject const & r) const
        {
            return r.cost_ < l.cost_;
        }
    };


    // provide copy constructor and assignment operator
    // to more realistically simulate a heterogeneous object 
    // where simple memset/memcopy is not possible
    HeapObject(const HeapObject &copy_from) 
	: super_type(copy_from), cost_(copy_from.cost_)
    {}

    inline HeapObject& operator= (const HeapObject &assign_from)
    {
	super_type::operator=(assign_from);
	cost_ = assign_from.cost_;
	return *this;
    }

};


template<int PAYLOADSIZE = 0, class CONVERSIONFUNCTOR=ConversionNOP>
struct SimulatedWatershedTests : TestData {
    enum { payload_size = PAYLOADSIZE };
    typedef size_t IdType; // FIXME: try vigra::MultiArrayIndex
    typedef size_t WeightType;
    typedef TestData::image_type image_type;

    // typedef IdType vertex_descriptor;
    typedef typename CONVERSIONFUNCTOR::vertex_descriptor vertex_descriptor;
    typedef typename CONVERSIONFUNCTOR::algo_vertex_descriptor algo_vertex_descriptor;
    typedef typename CONVERSIONFUNCTOR::queue_vertex_descriptor queue_vertex_descriptor;

    vector<vertex_descriptor> preparedDescriptorsFromLog;

    // constructor should prepare a cache of queue objects
    SimulatedWatershedTests() {
	// iterate log once, convert the ScanOrderIndices into the right objects
	cerr << "preparing log: ";
	bool type; WeightType priority; IdType id;
	//	preparedDescriptorsFromLog.resize(numentries);
	preparedDescriptorsFromLog.resize(numentries/2); // this assumes that id's are numberered 0...numentries/2-1 and each appears twice in log.
	TICTOCLOOP_BEGIN(1,1)
	reset();
	while (getNext(type, priority, id)) {
// 	    // FIXME: Only need to store data for pushes in fact!
// 	    preparedDescriptorsFromLog.push_back(CONVERSIONFUNCTOR::convert(id));
 	    preparedDescriptorsFromLog[id] = CONVERSIONFUNCTOR::convert(in, id);
	}	
	TICTOCLOOP_END
    }

#ifdef SIM_PQ_WATERSHED
    void simulatedWatershed() {
	cerr << "simulatedWatershed() using payload size "<< PAYLOADSIZE << " and conversion code from " << CONVERSIONFUNCTOR::name() << endl;
	//typedef HeapObject<WeightType, vertex_descriptor, payload_size> QueueObject;
	typedef HeapObject<WeightType, queue_vertex_descriptor, payload_size> QueueObject;
	algo_vertex_descriptor tmp_algo_vertex_descriptor;
	QueueObject tmp_queue_object(0, queue_vertex_descriptor());
	
	typedef std::priority_queue<QueueObject, std::vector<QueueObject>, typename QueueObject::Compare>  SeedRgNodeHeap;
	SeedRgNodeHeap pheap;
	bool type; WeightType priority; IdType id;
	TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	reset();
	while (getNext(type, priority, id)) {
	    if (type) { // push
		CONVERSIONFUNCTOR::execPush(in, preparedDescriptorsFromLog[id], tmp_queue_object.id_);
		tmp_queue_object.cost_ = priority;
		pheap.push(tmp_queue_object);
	    }
	    else { // pop
		const QueueObject &qo(pheap.top());
		if (0) {
		    // DEBUGGING ONLY: check if correct according to log
		    // (only works if same tiebreaker is used, of course.)
		    shouldEqual(qo.cost_, priority);
		    // shouldEqual(qo.id_, id);
		    // ==> Works OK (priorities OK) but: makes a huge difference in running time
		    // (@0: 483ms -> 6394ms)
		}
		CONVERSIONFUNCTOR::execPop(in, qo.id_, tmp_algo_vertex_descriptor);
		doNothing((void*)&tmp_algo_vertex_descriptor);
		pheap.pop();
	    }
	}
	TICTOCLOOP_END
	//cerr << counter << " queue operations processed." << endl;
	shouldEqual(counter, 4096766);
    }
#endif



#ifdef SIM_TURBO_WATERSHED
    void simulatedTurboWatershed() {
	cerr << "simulatedTurboWatershed() using payload size "<< PAYLOADSIZE << " and conversion code from " << CONVERSIONFUNCTOR::name() << endl;

	typedef TurboQueueObject<queue_vertex_descriptor, payload_size> QueueObject;	
	queue_vertex_descriptor tmp_queue_vertex_descriptor = queue_vertex_descriptor();
	QueueObject tmp_queue_object(tmp_queue_vertex_descriptor);
	algo_vertex_descriptor tmp_algo_vertex_descriptor;
	const size_t nGrayLevels = maxval + 1;  // maximum occuring in test dataset 

	// define one queue per gray level
	std::vector<std::queue<QueueObject> > queues(nGrayLevels);
	
	bool type; size_t priority; size_t id;
	TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	reset();
	while (getNext(type, priority, id)) {
	    if (type) { // push
		CONVERSIONFUNCTOR::execPush(in, preparedDescriptorsFromLog[id], tmp_queue_object.id_);
		queues[priority].push(tmp_queue_object);
	    }
	    else { // pop
		// should(queues[priority].size() > 0);
		const QueueObject &qo = queues[priority].front();
		// Do necessary conversions to iterate neighborhood HERE:
		CONVERSIONFUNCTOR::execPop(in, qo.id_, tmp_algo_vertex_descriptor);
		doNothing((void*)&tmp_algo_vertex_descriptor);
		queues[priority].pop();
	    }
	}
	TICTOCLOOP_END
	//cerr << counter << " queue operations processed." << endl;
	shouldEqual(counter, 4096766);
	// Check if all queues are empty by now
	for (size_t i=0; i < nGrayLevels; ++i)
	    shouldEqual(queues[i].size(), 0);
    }
#endif

};



struct GridGraphsTestSuite
: public vigra::test_suite
{
    GridGraphsTestSuite()
    : vigra::test_suite("GridGraphsTestSuite")
    {
#ifdef SIM_PQ_WATERSHED
	add(testCase((&SimulatedWatershedTests<0,ConversionNOP>::simulatedWatershed)));
// 	add(testCase(&SimulatedWatershedTests<2>::simulatedWatershed));
//  	add(testCase(&SimulatedWatershedTests<4>::simulatedWatershed));
//  	add(testCase(&SimulatedWatershedTests<8>::simulatedWatershed));
//	add(testCase((&SimulatedWatershedTests<16,ConversionNOP>::simulatedWatershed)));
//  	add(testCase(&SimulatedWatershedTests<32>::simulatedWatershed));
//  	add(testCase(&SimulatedWatershedTests<64>::simulatedWatershed));

	add(testCase((&SimulatedWatershedTests<0,SOI_OnQueue>::simulatedWatershed)));
	//	add(testCase((&SimulatedWatershedTests<16,SOI_OnQueue>::simulatedWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue>::simulatedWatershed)));
// 	add(testCase((&SimulatedWatershedTests<16,CoordsOnQueue>::simulatedWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue>::simulatedWatershed)));
//	add(testCase((&SimulatedWatershedTests<16,IterOnQueue>::simulatedWatershed)));
#endif // SIM_PQ_WATERSHED



#ifdef SIM_TURBO_WATERSHED
	// Turbo Watershed Simulation:
	add(testCase((&SimulatedWatershedTests<0,ConversionNOP>::simulatedTurboWatershed)));
// 	add(testCase(&SimulatedWatershedTests<2>::simulatedTurboWatershed));
//  	add(testCase(&SimulatedWatershedTests<4>::simulatedTurboWatershed));
//  	add(testCase(&SimulatedWatershedTests<8>::simulatedTurboWatershed));
//	add(testCase((&SimulatedWatershedTests<16,ConversionNOP>::simulatedTurboWatershed)));
//  	add(testCase(&SimulatedWatershedTests<32>::simulatedTurboWatershed));
//  	add(testCase(&SimulatedWatershedTests<64>::simulatedTurboWatershed));

	add(testCase((&SimulatedWatershedTests<0,SOI_OnQueue>::simulatedTurboWatershed)));
	//	add(testCase((&SimulatedWatershedTests<16,SOI_OnQueue>::simulatedTurboWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue>::simulatedTurboWatershed)));
// 	add(testCase((&SimulatedWatershedTests<16,CoordsOnQueue>::simulatedTurboWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue>::simulatedTurboWatershed)));
// 	add(testCase((&SimulatedWatershedTests<16,IterOnQueue>::simulatedTurboWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue_NoCopies>::simulatedTurboWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue_ManyNOPs>::simulatedTurboWatershed)));

#endif // SIM_TURBO_WATERSHED

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
