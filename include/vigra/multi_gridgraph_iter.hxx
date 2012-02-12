// version of the iterator that yields the vertex descriptor itself on
// dereference operations
// (used for BGL compatibility)

template <unsigned int N, class T>
class GridGraphBGLNodeIterator
: public StridedScanOrderIterator<N, T, T&, T*>
{
    typedef StridedScanOrderIterator<N, T, T&, T*> base_type;
    typedef GridGraphBGLNodeIterator<N, T> self_type;
  public:

    typedef typename base_type::value_type value_type;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::shape_type shape_type;
    typedef MultiArrayIndex difference_type;
    typedef GridGraphBGLNodeIterator iterator;
    typedef std::random_access_iterator_tag iterator_category;

    GridGraphBGLNodeIterator() : base_type()
    {}

    GridGraphBGLNodeIterator(const base_type &other)
    //	: base_type(const_cast<T*>(other.ptr()), other.shape(), other.strides()), base_type::index_(other.index())
	: base_type(other)
    {}

//     GridGraphBGLNodeIterator(self_type &other)
// 	: base_type(other.ptr(), other.shape(), other.strides())
//     {}
     GridGraphBGLNodeIterator(const self_type &other)
	 : base_type(other)//const_cast<T*>(other.ptr()), other.shape(), other.strides())
     {}

    GridGraphBGLNodeIterator(pointer i, 
			     shape_type const & shape, shape_type const & strides)
    : base_type(i, shape, strides)
    {}


    using base_type::point;
    using base_type::index;
    using base_type::neighborhoodType;
    using base_type::operator->;


    MultiArrayIndex operator*() const
    {
        return this->index_;
    }

    MultiArrayIndex operator[](MultiArrayIndex i)
    {
        self_type t(*this);
        t.moveToScanOrderIndex(this->index_+i);
        return *t;
    }

    MultiArrayIndex operator[](MultiArrayIndex i) const
    {
        self_type t(*this);
        t.moveToScanOrderIndex(this->index_+i);
        return *t;
    }


};

