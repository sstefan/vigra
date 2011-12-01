/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2011 by Stefan Schmidt and Ullrich Koethe         */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.3.0, Sep 10 2004 )                                    */
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


/*
 * multi_iterator_coupled.hxx
 *
 *  Created on: Nov 4, 2011
 *      Author: Stefan Schmidt
 */

#ifndef MULTI_ITERATOR_COUPLED_HXX_
#define MULTI_ITERATOR_COUPLED_HXX_

#include "multi_iterator.hxx"

namespace vigra {


struct NoImage {
    typedef void value_type;
    typedef void reference;
    typedef void pointer;
};



template <unsigned int N, class T, class REFERENCE, class POINTER>
class PerImageState {
public:
    typedef PerImageState<N, T, REFERENCE, POINTER> self_type;

    typedef T value_type;
    typedef POINTER pointer_type;
    typedef T const * const_pointer_type;
    typedef REFERENCE reference_type;
    typedef T const & const_reference_type;
    typedef typename MultiArrayShape<N>::type shape_type;

    PerImageState(POINTER pointer, const shape_type& strides)
	: pointer_(pointer), strides_(strides)
    {}

    template<unsigned int LEVEL>
    inline void increment() {
        pointer_ += strides_[LEVEL];
    }
    template<int LEVEL>
    inline void decrement() {
        pointer_ -= strides_[LEVEL];
    }
    // TODO: test if making the above a default case of the this hurts performance
    template<int LEVEL>
    inline void increment(const size_t multiple) {
        pointer_ += multiple*strides_[LEVEL];
    }
    template<int LEVEL>
    inline void decrement(const size_t multiple) {
        pointer_ -= multiple*strides_[LEVEL];
    }

    // ptr access
    reference_type operator*()
    {
        return *pointer_;
    }

    const_reference_type operator*() const
    {
        return *pointer_;
    }

    pointer_type operator->()
    {
        return pointer_;
    }

    const_pointer_type operator->() const
    {
        return pointer_;
    }

    pointer_type ptr()
    {
        return pointer_;
    }

    const_pointer_type ptr() const
    {
        return pointer_;
    }

    shape_type const & strides() const
    {
        return strides_;
    }

protected:
    pointer_type pointer_;
    shape_type strides_;
};


#if 1
// specialized template to handle no image (default) case
template<unsigned int N>
class PerImageState<N, NoImage, NoImage, NoImage> {
public:
	template<int LEVEL>
	inline void increment(const size_t multiple=1) {}
	template<int LEVEL>
	inline void decrement(const size_t multiple=1) {}
};
#endif




struct ChainEnd {
    typedef ChainEnd self_type;
    typedef ChainEnd next_type;

    template<int LEVEL>
    inline void increment(const size_t multiple=1) {}
    template<int LEVEL>
    inline void decrement(const size_t multiple=1) {}
};


template <unsigned int N, class T, class REFERENCE, class POINTER, class TAIL=ChainEnd>
class ChainedPerImageState
    : public PerImageState<N, T, REFERENCE, POINTER>
{
public:
    typedef PerImageState<N, T, REFERENCE, POINTER> base_type;
    typedef ChainedPerImageState<N, T, REFERENCE, POINTER, TAIL> self_type;
    typedef TAIL next_type;

    typedef T value_type;
    typedef POINTER pointer_type;
    typedef T const * const_pointer_type;
    typedef REFERENCE reference_type;
    typedef T const & const_reference_type;
    typedef typename MultiArrayShape<N>::type shape_type;

    ChainedPerImageState(POINTER pointer, const shape_type& strides, const next_type& next = ChainEnd())
	: base_type(pointer, strides),
	  next_(next)
    {}

    template<unsigned int LEVEL>
    inline void increment() {
	base_type::template increment<LEVEL>();
	next_.template increment<LEVEL>();
    }
    template<int LEVEL>
    inline void decrement() {
	base_type::template decrement<LEVEL>();
	next_.template decrement<LEVEL>();
    }
    // TODO: test if making the above a default case of the this hurts performance
    template<int LEVEL>
    inline void increment(const size_t multiple) {
	base_type::template increment<LEVEL>(multiple);
	next_.template increment<LEVEL>(multiple);
    }
    template<int LEVEL>
    inline void decrement(const size_t multiple) {
	base_type::template decrement<LEVEL>(multiple);
	next_.template decrement<LEVEL>(multiple);
    }

    next_type& next() 
    {
	return next_;
    }

protected:
    next_type next_;
};



namespace detail {

// Traits class to determine the type of a component in the chain.
// Also recursively selects the right component on query.
template<int B, class CHAIN, class COUPLEDITERATOR>
struct ComponentSelector {
    typedef typename ComponentSelector<B-1, typename CHAIN::next_type, COUPLEDITERATOR>::type type;

    static type& get(COUPLEDITERATOR &owner, CHAIN &chain) {
	return ComponentSelector<B-1, typename CHAIN::next_type, COUPLEDITERATOR>::get(owner, chain.next());
    }
};

#if 1 
template<class CHAIN, class COUPLEDITERATOR>
struct ComponentSelector<0, CHAIN, COUPLEDITERATOR> {
    typedef typename CHAIN::base_type type;

    static type& get(COUPLEDITERATOR &owner, CHAIN &chain) {
	return chain;
    }
};
#endif

template<int B, class COUPLEDITERATOR>
struct ComponentSelector<B, ChainEnd, COUPLEDITERATOR> {
    typedef COUPLEDITERATOR type;

    static type& get(COUPLEDITERATOR &owner, ChainEnd &chain) {
	return owner;
    }
};


#if 1
// FIXME: still necessary? -> Perhaps this type resolution could be solved with less template specializations.
template<class COUPLEDITERATOR>
struct ComponentSelector<0, ChainEnd, COUPLEDITERATOR> {
    typedef COUPLEDITERATOR type;

    static type& get(COUPLEDITERATOR &owner, ChainEnd &chain) {
	return owner;
    }
};
#endif


} // namespace detail






/********************************************************/
/*                                                      */
/*               CoupledScanOrderIterator<N>            */
/*                                                      */
/********************************************************/

/** \brief Iterate over multiple images simultaneously in scan order. Zero images is a special case;
    The coordinates can be accessed as a special band.

<b>\#include</b> \<vigra/multi_iterator.hxx\>

Namespace: vigra
*/

template <unsigned int N,
	  class CHAIN = PerImageState<N, NoImage, NoImage, NoImage>, // FIXME: Use chain as default!
	  unsigned int M = N>
class CoupledScanOrderIterator
#ifndef DOXYGEN  // doxygen doesn't understand this inheritance
: protected CoupledScanOrderIterator<N, CHAIN, M-1>
#endif
{
    typedef CoupledScanOrderIterator<N, CHAIN, M-1> base_type;
    enum { level = M-1 };

  public:

    typedef typename base_type::shape_type shape_type;
    typedef MultiArrayIndex difference_type;
    typedef CoupledScanOrderIterator iterator;
    typedef std::random_access_iterator_tag iterator_category;

    typedef typename base_type::PointerChain_type PointerChain_type;

    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference; // FIXME: do we need both?

    CoupledScanOrderIterator()
    {}

    CoupledScanOrderIterator(shape_type const & shape,
			     PointerChain_type pc = PointerChain_type())
    : base_type(shape, pc)
    {}

    CoupledScanOrderIterator & operator++()
    {
        base_type::operator++();
        if(this->point_[level-1] == this->shape_[level-1])
        {
            base_type::reset();

            this->pointers_.template increment<level>();

            ++this->point_[level];
        }
        return *this;
    }

    CoupledScanOrderIterator operator++(int)
    {
        CoupledScanOrderIterator res(*this);
        ++*this;
        return res;
    }

    CoupledScanOrderIterator & operator+=(MultiArrayIndex i)
    {
        this->moveToScanOrderIndex(this->index_+i);
        return *this;
    }

    CoupledScanOrderIterator & operator+=(const shape_type &coordOffset)
    {
        this->moveRelative(detail::CoordinateToScanOrder<N>::exec(this->shape_, coordOffset),
			   coordOffset);
        return *this;
    }

    CoupledScanOrderIterator & operator--()
    {
        base_type::operator--();
        if(this->point_[level-1] == -1)
        {
            base_type::inverseReset();
            this->pointers_.template decrement<level>();
            --this->point_[level];
        }
        return *this;
    }

    CoupledScanOrderIterator operator--(int)
    {
        CoupledScanOrderIterator res(*this);
        --*this;
        return res;
    }

    CoupledScanOrderIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    CoupledScanOrderIterator & operator-=(const shape_type &coordOffset)
    {
	return operator+=(-coordOffset);
    }

    CoupledScanOrderIterator getEndIterator() const
    {
        CoupledScanOrderIterator res(*this);
        res.moveToScanOrderIndex(prod(this->shape_));
        return res;
    }

    bool atBorder() const
    {
        return base_type::atBorder() ||
                this->point_[level] == 0 ||
                this->point_[level] == this->shape_[level] - 1;
    }

    unsigned int neighborhoodType() const
    {
        unsigned int res = base_type::neighborhoodType();
        if(this->point_[level] == 0)
            res |= (1 << 2*level);
        if(this->point_[level] == this->shape_[level]-1)
            res |= (2 << 2*level);
        return res;
    }

    CoupledScanOrderIterator operator+(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) += d;
    }

    CoupledScanOrderIterator operator-(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) -= d;
    }

    CoupledScanOrderIterator operator+(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) += coordOffset;
    }

    CoupledScanOrderIterator operator-(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) -= coordOffset;
    }


    MultiArrayIndex operator-(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator-(r);
    }

    bool operator==(CoupledScanOrderIterator const & r)
    {
        return base_type::operator==(r);
    }

    bool operator!=(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator!=(r);
    }

    bool operator<(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator<(r);
    }

    bool operator<=(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator<=(r);
    }

    bool operator>(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator>(r);
    }

    bool operator>=(CoupledScanOrderIterator const & r) const
    {
        return base_type::operator>=(r);
    }

    using base_type::point;
    using base_type::shape;
    using base_type::index;
    using base_type::operator[]; 
    using base_type::operator*;

    using base_type::top;
    using base_type::get;

  protected:
    void reset()
    {
        this->pointers_.template decrement<level>(this->shape_[level]);
        this->point_[level] = 0;
    }

    void inverseReset()
    {
        this->pointers_.template increment<level>(this->shape_[level]);
        this->point_[level] = this->shape_[level]-1;
    }


    CoupledScanOrderIterator & moveRelative(const MultiArrayIndex &indexOffset,
					    const shape_type &coordOffset)
    {
	base_type::moveRelative(indexOffset, coordOffset);
	this->pointers_.template increment<level>(coordOffset[level]);
	this->point_[level] += coordOffset[level];
        return *this;
    }
};



template <unsigned int N, class CHAIN>
class CoupledScanOrderIterator<N, CHAIN, 1>
{
    enum { level = 0 };

  public:

    typedef class CoupledScanOrderIterator<N, CHAIN, 1> this_type;

    typedef typename MultiArrayShape<N>::type shape_type;
    typedef MultiArrayIndex difference_type;
    typedef CoupledScanOrderIterator iterator;
    typedef std::random_access_iterator_tag iterator_category;

    typedef CHAIN PointerChain_type;

    // dereferencing this iterator itself yields the coordinate tuple:
    typedef shape_type const & reference;
    typedef shape_type const & const_reference; // FIXME: do we need both?
    

    CoupledScanOrderIterator()
	: 
	pointers_(PointerChain_type()),
	index_(0)
    {}

    CoupledScanOrderIterator(shape_type const & shape,
			     PointerChain_type pc = PointerChain_type())
    : shape_(shape),
      pointers_(pc),
      index_(0)
    {}

    CoupledScanOrderIterator & operator++()
    {
        pointers_.template increment<level>();
        ++point_[level];
        ++index_;
        return *this;
    }

    CoupledScanOrderIterator operator++(int)
    {
        CoupledScanOrderIterator res(*this);
        ++*this;
        return res;
    }

    CoupledScanOrderIterator & operator+=(MultiArrayIndex i)
    {
        this->moveToScanOrderIndex(index_+i);
        return *this;
    }

    CoupledScanOrderIterator  & operator+=(const shape_type &coordOffset)
    {
        moveRelative(detail::CoordinateToScanOrder<N>::exec(shape_, coordOffset),
		     coordOffset);
        return *this;
    }
    CoupledScanOrderIterator & operator-=(const shape_type &coordOffset)
    {
	return operator+=(-coordOffset);
    }

    CoupledScanOrderIterator & operator--()
    {
        pointers_.template decrement<level>();
        --point_[level];
        --index_;
        return *this;
    }

    CoupledScanOrderIterator operator--(int)
    {
        CoupledScanOrderIterator res(*this);
        --this;
        return res;
    }

    CoupledScanOrderIterator & operator-=(MultiArrayIndex i)
    {
        return operator+=(-i);
    }

    reference operator[](MultiArrayIndex i)
    {
        CoupledScanOrderIterator t(*this);
        t.moveToScanOrderIndex(index_+i);
        return *t;
    }

    const_reference operator[](MultiArrayIndex i) const
    {
        CoupledScanOrderIterator t(*this);
        t.moveToScanOrderIndex(index_+i);
        return *t;
    }

    CoupledScanOrderIterator
    operator+(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) += d;
    }

    CoupledScanOrderIterator
    operator-(MultiArrayIndex d) const
    {
        return CoupledScanOrderIterator(*this) -= d;
    }

    CoupledScanOrderIterator operator+(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) += coordOffset;
    }
    CoupledScanOrderIterator operator-(const shape_type &coordOffset) const
    {
        return CoupledScanOrderIterator(*this) -= coordOffset;
    }

    MultiArrayIndex
    operator-(CoupledScanOrderIterator const & r) const
    {
        return index() - r.index();
    }

    bool
    operator==(CoupledScanOrderIterator const & r)
    {
        return index() == r.index();
    }

    bool
    operator!=(CoupledScanOrderIterator const & r) const
    {
        return index() != r.index();
    }

    bool
    operator<(CoupledScanOrderIterator const & r) const
    {
        return index() < r.index();
    }

    bool
    operator<=(CoupledScanOrderIterator const & r) const
    {
        return index() <= r.index();
    }

    bool
    operator>(CoupledScanOrderIterator const & r) const
    {
        return index() > r.index();
    }

    bool
    operator>=(CoupledScanOrderIterator const & r) const
    {
        return index() >= r.index();
    }

    bool atBorder() const
    {
        return point_[level] == 0 || point_[level] == shape_[level] - 1;
    }

    MultiArrayIndex index() const
    {
        return index_;
    }

    shape_type const & point() const
    {
        return point_;
    }

    //! operator* overridden, provides coordinates
    shape_type const & operator*()
    {
        return point_;
    }

    shape_type const & shape() const
    {
        return shape_;
    }

    CoupledScanOrderIterator getEndIterator() const
    {
        CoupledScanOrderIterator res(*this);
        res.moveToScanOrderIndex(prod(shape_));
        return res;
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


    PointerChain_type& top() { return pointers_; }

    template<unsigned int B> 
    typename detail::ComponentSelector<B, PointerChain_type, this_type>::type&
    get() 
    {
	return detail::ComponentSelector<B, PointerChain_type, this_type>::get(*this, pointers_);
    }

  protected:
    void reset()
    {
        this->pointers_.template decrement<level>(shape_[level]);
        point_[level] = 0;
    }

    void inverseReset()
    {
        this->pointers_.template increment<level>(shape_[level]);
        point_[level] = shape_[level] - 1;
    }

    void moveToScanOrderIndex(MultiArrayIndex newIndex)
    {
        index_ = newIndex;
        detail::MoveToScanOrderIndex<N-1>::exec(newIndex, shape_, point_, pointers_);
    }

    CoupledScanOrderIterator & moveRelative(const MultiArrayIndex &indexOffset,
					    const shape_type &coordOffset)
    {
	point_[level] += coordOffset[level];
	
        index_+= indexOffset;
	pointers_.template increment<level>(coordOffset[level]);
	
        return *this;
    }   
    
    shape_type point_, shape_; // FIXME: shape_ should become upperleft_, lowerright_ tuple
    PointerChain_type pointers_;
    MultiArrayIndex index_;
};









// --------------------------------------------------
// Factory



// convenience adapter for MultiArrayViews:
template <class VIEW>
struct PtrHolder {
    typedef PerImageState<VIEW::actual_dimension,
			  typename VIEW::value_type,
			  typename VIEW::reference,
			  typename VIEW::pointer> type;
    typedef PerImageState<VIEW::actual_dimension,
			  typename VIEW::value_type,
			  typename VIEW::const_reference,
			  typename VIEW::const_pointer> const_type;
};


// convenience adapter for MultiArrayViews:
template <class VIEW, class TAIL=ChainEnd>
struct ChainedPtrHolder {
    typedef ChainedPerImageState<VIEW::actual_dimension,
				 typename VIEW::value_type,
				 typename VIEW::reference,
				 typename VIEW::pointer,
				 TAIL> type;
};
template <class TAIL>
struct ChainedPtrHolder<NoImage, TAIL> {
    typedef ChainEnd type;
};




template<int N,
	 class I1=NoImage,
	 class I2=NoImage,
	 class I3=NoImage,
	 class I4=NoImage,
	 class I5=NoImage>
class CoupledScanOrderIteratorFactory 
{
public:
    typedef typename MultiArrayShape<N>::type shape_type;
    typedef typename ChainedPtrHolder<I1, 
	    typename ChainedPtrHolder<I2, 
	    typename ChainedPtrHolder<I3, 
	    typename ChainedPtrHolder<I4, 
	    typename ChainedPtrHolder<I5>::type >::type >::type >::type >::type Chain;

    typedef CoupledScanOrderIterator<N, Chain> 
        coupled_iter_type;

    static coupled_iter_type makeCoupledIterator(const shape_type &shape) {
	return coupled_iter_type(shape);
    }
    static coupled_iter_type makeCoupledIterator(I1 &i1) {
	return coupled_iter_type(i1.shape(), 
				 Chain(i1.data(), i1.stride()));
    }
    static coupled_iter_type makeCoupledIterator(I1 &i1, I2 &i2) {
	return coupled_iter_type(i1.shape(), 
				 Chain(i1.data(), i1.stride(), 
			       typename Chain::next_type(i2.data(), i2.stride())));
    }
    static coupled_iter_type makeCoupledIterator(I1 &i1, I2 &i2, I3 &i3) {
	return coupled_iter_type(i1.shape(), 
				 Chain(i1.data(), i1.stride(), 
				       typename Chain::next_type(i2.data(), i2.stride(),
							 typename Chain::next_type::next_type(i3.data(), i3.stride()))));
    }
    static coupled_iter_type makeCoupledIterator(I1 &i1, I2 &i2, I3 &i3, I4 &i4) {
	return coupled_iter_type(i1.shape(), 
	    Chain(i1.data(), i1.stride(), 
	    typename Chain::next_type(i2.data(), i2.stride(),
	    typename Chain::next_type::next_type(i3.data(), i3.stride(),
	    typename Chain::next_type::next_type::next_type(i4.data(), i4.stride())))));
    }
    static coupled_iter_type makeCoupledIterator(I1 &i1, I2 &i2, I3 &i3, I4 &i4, I5 &i5) {
	return coupled_iter_type(i1.shape(), 
	    Chain(i1.data(), i1.stride(), 
	    typename Chain::next_type(i2.data(), i2.stride(),
	    typename Chain::next_type::next_type(i3.data(), i3.stride(),
	    typename Chain::next_type::next_type::next_type(i4.data(), i4.stride(),
	    typename Chain::next_type::next_type::next_type::next_type(i5.data(), i5.stride()))))));
    }

};


}; // namespace vigra

#endif /* MULTI_ITERATOR_COUPLED_HXX_ */
