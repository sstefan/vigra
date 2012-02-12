template <unsigned int N, class T, class Flags = UInt8, unsigned int M = N>
class GridGraphNodeIterator
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: protected GridGraphNodeIterator<N, T, Flags, M-1>
#endif
{
    typedef GridGraphNodeIterator<N, T, Flags, M-1> base_type;
  public:
    enum { level = M-1 };

    typedef typename base_type::value_type value_type;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::flags_pointer flags_pointer;
    typedef typename base_type::shape_type shape_type;

    GridGraphNodeIterator(pointer i, flags_pointer f,
                          shape_type const & shape, shape_type const & strides)
    : base_type(i, f, shape, strides)
    {}

    GridGraphNodeIterator & operator++()
    {
        base_type::operator++();
        if(this->point_[level-1] == this->shape_[level-1] && 
           this->point_[level] < this->shape_[level])
        {
            base_type::reset();
            this->i_ += this->strides_[level];
            ++this->point_[level];
        }
        return *this;
    }

        /** Advance to next starting location.
         */
    GridGraphNodeIterator operator++(int)
    {
        GridGraphNodeIterator res(*this);
        ++*this;
        return res;
    }

    GridGraphNodeIterator & operator+=(MultiArrayIndex i)
    {
        this->moveToId(this->index_+i);
        return *this;
    }

    GridGraphNodeIterator & operator--()
    {
        base_type::operator--();
        if(this->point_[level-1] == -1 && this->point_[level] > 0)
        {
            base_type::inverseReset();
            this->i_ -= this->strides_[level];
            --this->point_[level];
        }
        return *this;
    }

        /** Advance to next starting location.
         */
    GridGraphNodeIterator operator--(int)
    {
        GridGraphNodeIterator res(*this);
        --*this;
        return res;
    }

    GridGraphNodeIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    GridGraphNodeIterator getEndIterator() const
    {
        GridGraphNodeIterator res(*this);
        res.makeEndIterator();
        return res;
    }
    
    using base_type::point;
    using base_type::index;
    using base_type::atBorder;
    using base_type::operator*;
    using base_type::operator->;
    using base_type::operator[];
    using base_type::operator==;
    using base_type::operator!=;
    
    unsigned int neighborhoodType() const
    {
        unsigned int res = base_type::neighborhoodType();
        if(this->point_[level] == 0)
            res |= (1 << 2*level);
        if(this->point_[level] == this->shape_[level]-1)
            res |= (2 << 2*level);
        return res;
    }

  protected:
    void reset()
    {
        this->i_ -= this->shape_[level]*this->strides_[level];
        this->point_[level] = 0;
    }

    void inverseReset()
    {
        this->i_ += this->shape_[level]*this->strides_[level];
        this->point_[level] = this->shape_[level] - 1;
    }    
};

template <unsigned int N, class T, class Flags>
class GridGraphNodeIterator<N, T, Flags, 1>
{
  public:
    enum { level = 0 };

    typedef T value_type;
    typedef T * pointer;
    typedef T const * const_pointer;
    typedef T & reference;
    typedef T const & const_reference;
    typedef Flags * flags_pointer;
    typedef typename MultiArrayShape<N>::type shape_type;

    GridGraphNodeIterator(pointer i, flags_pointer f,
                          shape_type const & shape, shape_type const & strides)
    : i_(i),
      f_(f),
      shape_(shape),
      strides_(strides),
      index_(0)
    {}

    GridGraphNodeIterator & operator++()
    {
        i_ += strides_[level];
        ++point_[level];
        ++index_;
        return *this;
    }

    GridGraphNodeIterator operator++(int)
    {
        GridGraphNodeIterator res(*this);
        ++*this;
        return res;
    }

    GridGraphNodeIterator & operator+=(MultiArrayIndex i)
    {
        this->moveToId(index_+i);
        return *this;
    }

    GridGraphNodeIterator & operator--()
    {
        i_ -= strides_[level];
        --point_[level];
        --index_;
        return *this;
    }

    GridGraphNodeIterator operator--(int)
    {
        GridGraphNodeIterator res(*this);
        --*this;
        return res;
    }

    GridGraphNodeIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }
    
    reference operator*()
    {
        return *i_;
    }

    const_reference operator*() const
    {
        return *i_;
    }

    pointer operator->()
    {
        return i_;
    }

    const_pointer operator->() const
    {
        return i_;
    }

    reference operator[](MultiArrayIndex i)
    {
        GridGraphNodeIterator t(*this);
        t.moveToId(index_+i);
        return *t;
    }

    const_reference operator[](MultiArrayIndex i) const
    {
        GridGraphNodeIterator t(*this);
        t.moveToId(index_+i);
        return *t;
    }

    bool atBorder() const
    {
        return f_[index_] != 0;
    }
    
    MultiArrayIndex index() const
    {
        return index_;
    }
    
    shape_type const & point() const
    {
        return point_;
    }
    
    unsigned int neighborhoodType() const
    {
        unsigned int res = 0;
        if(this->point_[level] == 0)
            res |= 1;
        if(this->point_[level] == this->shape_[level]-1)
            res |= 2;
        return res;
    }

    template <unsigned int K>
    bool operator==(GridGraphNodeIterator<N, T, Flags, K> const & other) const
    {
        return index_ == other.index();
    }
    
    template <unsigned int K>
    bool operator!=(GridGraphNodeIterator<N, T, Flags, K> const & other) const
    {
        return index_ != other.index();
    }
     
    GridGraphNodeIterator getEndIterator() const
    {
        GridGraphNodeIterator res(*this);
        res.makeEndIterator();
        return res;
    }
   
  protected:
    void reset()
    {
        i_ -= shape_[level]*strides_[level];
        point_[level] = 0;
    }
    
    void inverseReset()
    {
        i_ += shape_[level]*strides_[level];
        point_[level] = shape_[level] - 1;
    }
    
    void makeEndIterator()
    {
        point_ = shape_;
        index_ = prod(shape_);
    }

    void moveToId(MultiArrayIndex newId)
    {
        index_ = newId;
        detail::MoveToScanOrderIndex<N-1>::exec(newId, 1, i_, point_, shape_, strides_);
    }
    
    pointer i_;
    flags_pointer f_;
    shape_type point_, shape_, strides_;
    MultiArrayIndex index_;
};

