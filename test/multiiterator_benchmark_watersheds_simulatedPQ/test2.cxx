#define SIM_PQ_WATERSHED
#define SIM_TURBO_WATERSHED
#define SIM_PQ_WATERSHED_WITH_ALLOCATOR
// #define USE_COMPILED_TEST_DATA

// #define PAYLOAD_COPY_AND_INIT_LOOPS

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
typedef IMG::view_type IMGVIEW;


namespace std {
ostream& operator<<(ostream &stream, const IMGVIEW::iterator iter) {
    stream << iter.index() << "," << iter.point();
    return stream;
}
};

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
    typedef size_t queue_vertex_descriptor;
    typedef size_t algo_vertex_descriptor;

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	doNothing((void*)&in);
	doNothing((void*)&out); // redundand: done outside anyways.
    }

    static inline const queue_vertex_descriptor& execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	doNothing((void*)&in);
	doNothing((void*)&out); // pushed on queue. therefore not necessary
	return in;
    }

    static inline algo_vertex_descriptor convert(const IMG &img, const size_t &id)
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
    typedef vigra::MultiArrayIndex queue_vertex_descriptor;

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	doNothing((void*)&in);
	// on Pop, need to reconstruct coords
	// or even full SSOI to iterate neighborhood:
	out = static_cast<view_type>(img).begin() + in;
	doNothing((void*)&out); // redundand: done outside anyways.
    }


    static inline const queue_vertex_descriptor& execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	doNothing((void*)&in);
	// on Push, just extract SOI
	out = in.index();
	doNothing((void*)&out); // pushed on queue. therefore not necessary
	return out;
    }


    static inline algo_vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	algo_vertex_descriptor tmp = static_cast<view_type>(img).begin();
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

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    { 
	doNothing((void*)&in);
	// convert coords to full scan order iterator:
	out = static_cast<view_type>(img).begin() + img.coordinateToScanOrderIndex(in);
	doNothing((void*)&out); // redundand: done outside anyways.
    }

    static inline const queue_vertex_descriptor& execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    { 	// push executes for every unlabeled neighbor of a point.
	// we just store the (already available) coordinates here.
	doNothing((void*)&in);
	out = in.point();
	doNothing((void*)&out); // pushed on queue. therefore not necessary
	return out;
    }

    static inline algo_vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	algo_vertex_descriptor tmp = static_cast<view_type>(img).begin();
	tmp += id;
	return tmp;
    }

    static const string name() { return "CoordsOnQueue"; }
};




//! This class simulates the case that only the coordinates tuple is put on the queues,
//  and a full scanorderiterator is used in the algorithm:
struct CoordsOnQueue_OpPlus : ConversionNOP {
    typedef IMG::view_type view_type;

    // uses the coordinates on the queue:
    typedef IMG::difference_type  queue_vertex_descriptor;
    // algorithm uses the full scridedscanorderiterator:
    typedef view_type::iterator algo_vertex_descriptor; // the StridedScanOrderIterator

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    { 
	doNothing((void*)&in);
	// convert coords to full scan order iterator:
	//	out = static_cast<view_type>(img).begin() + img.coordinateToScanOrderIndex(in);
 	out = static_cast<view_type>(img).begin(); 
 	out += in;
	doNothing((void*)&out); // redundand: done outside anyways.
    }

    static inline const queue_vertex_descriptor& execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    { 	// push executes for every unlabeled neighbor of a point.
	// we just store the (already available) coordinates here.
	doNothing((void*)&in);
	out = in.point();
	doNothing((void*)&out); // pushed on queue. therefore not necessary
	return out;
    }

    static inline algo_vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	algo_vertex_descriptor tmp = static_cast<view_type>(img).begin();
	tmp += id;

	if (0) {
	    // some tests -> moved to test/multiarray/test.cxx
	    algo_vertex_descriptor tmp2 = static_cast<view_type>(img).begin();
	    tmp2 += tmp.point();
	    
	    shouldEqual(tmp2, tmp);
	
	    tmp2 -= tmp.point();
	    shouldEqual(tmp2, static_cast<view_type>(img).begin());
	    
	    tmp2 = static_cast<view_type>(img).begin() + tmp.point();
	    shouldEqual(tmp2, tmp);

	    tmp2 = tmp - tmp.point();
	    shouldEqual(tmp2, static_cast<view_type>(img).begin());
	}

	return tmp;
    }

    static const string name() { return "CoordsOnQueue_OpPlus"; }
};









//! This class simulates the case that a full StridedScanOrderIterator
//  is put on the queue.
//   In addition to the coordinates tuple, it already has the 
//   neighborhood information and can iterate its neighbors,
//   and furthermore the linear index / ptr is also available directly.
struct IterOnQueue : ConversionNOP {
    typedef IMG::view_type view_type;
    typedef view_type::iterator queue_vertex_descriptor;
    typedef view_type::iterator algo_vertex_descriptor; 

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	doNothing((void*)&in);
	// no special actions to take on pop. 
	//  however, we need to copy the object off the heap/queue:
	out = in;
	doNothing((void*)&out); // redundand: done outside anyways.
    }

    static inline const queue_vertex_descriptor& execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	doNothing((void*)&in);
	out = in; // FIXME: Can we avoid this copy/assignment as well?
	// FIXME: Instead of using the output argument, may directly pass reference to
	// input as return value. See below.
	doNothing((void*)&out); // pushed on queue. therefore not necessary
	return out;
    }

    static inline algo_vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	algo_vertex_descriptor res = static_cast<view_type>(img).begin();
	res += id; 
	return res;
    }

    static const string name() { return "IterOnQueue"; }
};


//! This class simulates the case that a full StridedScanOrderIterator
//  is put on the queue.
//  This variant should exactly cost as much as a NOP conversion with a dummy payload of corresponding size.
struct IterOnQueue_NoCopies : IterOnQueue {
    typedef IMG::view_type view_type;
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

    static inline const queue_vertex_descriptor&  execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	doNothing((void*)&in);	// do not even copy.
	doNothing((void*)&out);	// do not even copy.
	return in;
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

    static inline const queue_vertex_descriptor& execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
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
	return in;
    }

    static const string name() { return "IterOnQueue_ManyNOPs"; }
};









//! This class simulates the case that a reduced StridedScanOrderIterator 
//   (only its pointer, index and coords, but not shape and strides)
//   is put on the queue, as that information is redundant.
struct SlicedIterOnQueue : ConversionNOP {
    typedef IMG::view_type view_type;
    
    template<class SSOI>
    struct PartialStateOnQueue {
	typedef PartialStateOnQueue<SSOI> this_type;

	typename SSOI::pointer i_;
	typename SSOI::shape_type point_;
	MultiArrayIndex index_;

	PartialStateOnQueue() {}

	PartialStateOnQueue(SSOI& in) :
	    i_(in.ptr()), point_(in.point()), index_(in.index())
	{}


	struct MySSOI : public SSOI {
	    MySSOI(const SSOI& tmpl, const this_type &state)
		: SSOI(tmpl)
	    {
		this->i_ = state.i_;
		this->point_ = state.point_;
		this->index_ = state.index_;
	    }
	};

	SSOI toSSOI(const SSOI& tmpl) const {
	    MySSOI ret(tmpl, *this);
	    return ret;
	}
    };

    typedef view_type::iterator SSOI;
    typedef PartialStateOnQueue<SSOI> queue_vertex_descriptor;
    typedef view_type::iterator algo_vertex_descriptor; 

    static inline void execPop(const IMG &img, const queue_vertex_descriptor &in, algo_vertex_descriptor &out) 
    {
	// on pop, need to reconstruct a full descriptor from the state information:
	out = in.toSSOI(out);
    }

    static inline const queue_vertex_descriptor& execPush(const IMG &img, algo_vertex_descriptor &in, queue_vertex_descriptor &out)
    {
	out = queue_vertex_descriptor(in);
	return out;
    }

    static inline algo_vertex_descriptor convert(const IMG & img, const size_t &id)
    {
	algo_vertex_descriptor res = static_cast<view_type>(img).begin();
	res += id; 
	return res;
    }

    static const string name() { return "SlicedIterOnQueue"; }
};



template<class IDTYPE, int PAYLOADSIZE>
struct TurboQueueObject 
{
    IDTYPE id_;
    unsigned char data_[PAYLOADSIZE];
    TurboQueueObject(const IDTYPE &id)
	: id_(id) {
#ifdef PAYLOAD_COPY_AND_INIT_LOOPS
	for (size_t i=0; i<PAYLOADSIZE; i++)
	    data_[i] = 0;
#endif
    }

#ifdef PAYLOAD_COPY_AND_INIT_LOOPS
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
#endif
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
        bool operator()(HeapObject const * l,
                        HeapObject const * r) const
        {
            return r->cost_ < l->cost_;
        }
    };
    struct FakeCompareFast
    {
        // must implement > since priority_queue looks for largest element
        bool operator()(HeapObject const & l,
                        HeapObject const & r) const
        {
	    return false; // suppresses queue sorting effect
        }
        bool operator()(HeapObject const * l,
                        HeapObject const * r) const
        {
	    return false;
        }
    };
    struct FakeCompareSlow
    {
        // must implement > since priority_queue looks for largest element
        bool operator()(HeapObject const & l,
                        HeapObject const & r) const
        {
	    return true;  // reordering required all the time
        }
        bool operator()(HeapObject const * l,
                        HeapObject const * r) const
        {
	    return true;
        }
    };


#ifdef PAYLOAD_COPY_AND_INIT_LOOPS
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
#endif
};


// --------------------------------------------------
// Allocator 
template<class OBJ>
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

    OBJ * create(const OBJ& copyfrom)
    {
	if(!freelist_.empty())
            {
                OBJ * res = freelist_.top();
                freelist_.pop();
                *res = copyfrom; 
                return res;
            }
	
	return new OBJ(copyfrom);
    }

    void dismiss(OBJ * p)
    {
	freelist_.push(p);
    }

    std::stack<OBJ *> freelist_; 
};


// /Allocator 
// --------------------------------------------------



template<int PAYLOADSIZE = 0, class CONVERSIONFUNCTOR=ConversionNOP>
struct SimulatedWatershedTests : TestData {
    enum { payload_size = PAYLOADSIZE };
    typedef size_t IdType; // FIXME: try vigra::MultiArrayIndex
    typedef size_t WeightType;
    typedef TestData::image_type image_type;

    typedef typename CONVERSIONFUNCTOR::algo_vertex_descriptor vertex_descriptor;
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
	typedef HeapObject<WeightType, queue_vertex_descriptor, payload_size> QueueObject;
	cerr << "simulatedWatershed() using payload size "<< PAYLOADSIZE 
	     << " (QO size=" << sizeof(QueueObject) << ")"
	     << " and conversion code from " << CONVERSIONFUNCTOR::name() << endl;
	//typedef HeapObject<WeightType, vertex_descriptor, payload_size> QueueObject;
	algo_vertex_descriptor tmp_algo_vertex_descriptor;
	QueueObject tmp_queue_object(0, queue_vertex_descriptor());

	if (1)
	{ // the real measurement:
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
	if (1)
	{ // run with faked comparison operator (always true/always false) to single out effect of queue resorting
	    cerr << "    worst case queue sorting: ";
	    typedef std::priority_queue<QueueObject, std::vector<QueueObject>, typename QueueObject::FakeCompareSlow>  SeedRgNodeHeap;
	    // typedef std::priority_queue<QueueObject, std::vector<QueueObject>, typename QueueObject::FakeCompareFast>  SeedRgNodeHeap;
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
		shouldEqual(counter, 4096766);
	}
	if (1)
	{ // run with faked comparison operator (always true/always false) to single out effect of queue resorting
	    cerr << "    best case queue sorting:  ";
	    typedef std::priority_queue<QueueObject, std::vector<QueueObject>, typename QueueObject::FakeCompareFast>  SeedRgNodeHeap;
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
		shouldEqual(counter, 4096766);
	}
	    

    }
#endif


#ifdef SIM_PQ_WATERSHED_WITH_ALLOCATOR
    void simulatedWatershedAllocator() {
	typedef HeapObject<WeightType, queue_vertex_descriptor, payload_size> QueueObject;
	cerr << "simulatedWatershedALLOCATOR() using payload size "<< PAYLOADSIZE 
	     << " (QO size=" << sizeof(QueueObject) << ")"
	     << " and conversion code from " << CONVERSIONFUNCTOR::name() << endl;
	//typedef HeapObject<WeightType, vertex_descriptor, payload_size> QueueObject;
	algo_vertex_descriptor tmp_algo_vertex_descriptor;
	QueueObject tmp_queue_object(0, queue_vertex_descriptor());
	
	Allocator<QueueObject> allocator;

	typedef std::priority_queue<QueueObject*, std::vector<QueueObject *>, typename QueueObject::Compare>  SeedRgNodeHeap;
	SeedRgNodeHeap pheap;
	bool type; WeightType priority; IdType id;
	TICTOCLOOP_BEGIN(repetitions,best_of_repetitions)
	reset();
	while (getNext(type, priority, id)) {
	    if (type) { // push
		CONVERSIONFUNCTOR::execPush(in, preparedDescriptorsFromLog[id], tmp_queue_object.id_);
		tmp_queue_object.cost_ = priority;		
		pheap.push(allocator.create(tmp_queue_object));
	    }
	    else { // pop
		const QueueObject &qo(*pheap.top());
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
		allocator.dismiss(pheap.top());
		pheap.pop();
	    }
	}
	// free temporary memory
	while(pheap.size() != 0) {
	    allocator.dismiss(pheap.top());
	    pheap.pop();
	}
	TICTOCLOOP_END
	//cerr << counter << " queue operations processed." << endl;
	shouldEqual(counter, 4096766);
    }
#endif



#ifdef SIM_TURBO_WATERSHED
    void simulatedTurboWatershed() {
	typedef TurboQueueObject<queue_vertex_descriptor, payload_size> QueueObject;	
	cerr << "simulatedTURBOWatershed() using payload size "<< PAYLOADSIZE 
	     << " (QO size=" << sizeof(QueueObject) << ")"
	     << " and conversion code from " << CONVERSIONFUNCTOR::name() << endl;

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
#if 0
		CONVERSIONFUNCTOR::execPush(in, preparedDescriptorsFromLog[id], tmp_queue_object.id_);
		queues[priority].push(tmp_queue_object);
#else
		queues[priority].push(CONVERSIONFUNCTOR::execPush(in, preparedDescriptorsFromLog[id], tmp_queue_object.id_));
#endif
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
	add(testCase((&SimulatedWatershedTests<16,ConversionNOP>::simulatedWatershed)));
  	add(testCase(&SimulatedWatershedTests<32>::simulatedWatershed));
  	add(testCase(&SimulatedWatershedTests<80>::simulatedWatershed));

	add(testCase((&SimulatedWatershedTests<0,SOI_OnQueue>::simulatedWatershed)));
	//	add(testCase((&SimulatedWatershedTests<16,SOI_OnQueue>::simulatedWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue>::simulatedWatershed)));
 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue_OpPlus>::simulatedWatershed)));
// 	add(testCase((&SimulatedWatershedTests<16,CoordsOnQueue>::simulatedWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue>::simulatedWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue_ManyNOPs>::simulatedWatershed)));
//	add(testCase((&SimulatedWatershedTests<16,IterOnQueue>::simulatedWatershed)));
#endif // SIM_PQ_WATERSHED

#ifdef SIM_PQ_WATERSHED_WITH_ALLOCATOR
	add(testCase((&SimulatedWatershedTests<0,ConversionNOP>::simulatedWatershedAllocator)));
  	add(testCase((&SimulatedWatershedTests<80,ConversionNOP>::simulatedWatershedAllocator)));

	add(testCase((&SimulatedWatershedTests<0,SOI_OnQueue>::simulatedWatershedAllocator)));
	//	add(testCase((&SimulatedWatershedTests<16,SOI_OnQueue>::simulatedWatershedAllocator)));

 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue>::simulatedWatershedAllocator)));
 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue_OpPlus>::simulatedWatershedAllocator)));
// 	add(testCase((&SimulatedWatershedTests<16,CoordsOnQueue>::simulatedWatershedAllocator)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue>::simulatedWatershedAllocator)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue_ManyNOPs>::simulatedWatershedAllocator)));
//	add(testCase((&SimulatedWatershedTests<16,IterOnQueue>::simulatedWatershedAllocator)));
#endif // SIM_PQ_WATERSHED_WITH_ALLOCATOR


#ifdef SIM_TURBO_WATERSHED
	// Turbo Watershed Simulation:
	add(testCase((&SimulatedWatershedTests<0,ConversionNOP>::simulatedTurboWatershed)));
	add(testCase((&SimulatedWatershedTests<16,ConversionNOP>::simulatedTurboWatershed)));
	add(testCase((&SimulatedWatershedTests<32,ConversionNOP>::simulatedTurboWatershed)));
	add(testCase((&SimulatedWatershedTests<80,ConversionNOP>::simulatedTurboWatershed)));
// 	add(testCase(&SimulatedWatershedTests<2>::simulatedTurboWatershed));
//  	add(testCase(&SimulatedWatershedTests<4>::simulatedTurboWatershed));
//  	add(testCase(&SimulatedWatershedTests<8>::simulatedTurboWatershed));
//	add(testCase((&SimulatedWatershedTests<16,ConversionNOP>::simulatedTurboWatershed)));
//  	add(testCase(&SimulatedWatershedTests<32>::simulatedTurboWatershed));
//  	add(testCase(&SimulatedWatershedTests<64>::simulatedTurboWatershed));

	add(testCase((&SimulatedWatershedTests<0,SOI_OnQueue>::simulatedTurboWatershed)));
	//	add(testCase((&SimulatedWatershedTests<16,SOI_OnQueue>::simulatedTurboWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue>::simulatedTurboWatershed)));
 	add(testCase((&SimulatedWatershedTests<0,CoordsOnQueue_OpPlus>::simulatedTurboWatershed)));
// 	add(testCase((&SimulatedWatershedTests<16,CoordsOnQueue>::simulatedTurboWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,IterOnQueue>::simulatedTurboWatershed)));
// 	add(testCase((&SimulatedWatershedTests<16,IterOnQueue>::simulatedTurboWatershed)));

 	add(testCase((&SimulatedWatershedTests<0,SlicedIterOnQueue>::simulatedTurboWatershed)));

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
