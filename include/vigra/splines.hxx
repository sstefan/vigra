/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_SPLINES_HXX
#define VIGRA_SPLINES_HXX

#include <cmath>
#include "config.hxx"
#include "mathutil.hxx"
#include "polynomial.hxx"
#include "array_vector.hxx"
#include "fixedpoint.hxx"
#include "multi_array.hxx"
#include "linear_solve.hxx"
#include <stdexcept>

namespace vigra {

namespace autodiff {

template <class T, int N>
class DualVector;

} // namespace autodiff

/** \addtogroup MathFunctions Mathematical Functions
*/
//@{
/* B-Splines of arbitrary order and interpolating Catmull/Rom splines.

    <b>\#include</b> \<vigra/splines.hxx\><br>
    Namespace: vigra
*/
#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

/** Basic interface of the spline functors.

    Implements the spline functions defined by the recursion

    \f[ B_0(x) = \left\{ \begin{array}{ll}
                                  1 & -\frac{1}{2} \leq x < \frac{1}{2} \\
                                  0 & \mbox{otherwise}
                        \end{array}\right.
    \f]

    and

    \f[ B_n(x) = B_0(x) * B_{n-1}(x)
    \f]

    where * denotes convolution, and <i>n</i> is the spline order given by the
    template parameter <tt>ORDER</tt>. These spline classes can be used as
    unary and binary functors, as kernels for \ref resamplingConvolveImage(),
    and as arguments for \ref vigra::SplineImageView. Note that the spline order
    is given as a template argument.

    <b>\#include</b> \<vigra/splines.hxx\><br>
    Namespace: vigra
*/
template <int ORDER, class T = double>
class BSplineBase
{
  public:

        /** the value type if used as a kernel in \ref resamplingConvolveImage().
        */
    typedef T            value_type;
        /** the functor's unary argument type
        */
    typedef T            argument_type;
        /** the functor's first binary argument type
        */
    typedef T            first_argument_type;
        /** the functor's second binary argument type
        */
    typedef unsigned int second_argument_type;
        /** the functor's result type (unary and binary)
        */
    typedef T            result_type;
        /** the spline order
        */
    enum StaticOrder { order = ORDER };

        /** Create functor for given derivative of the spline. The spline's order
            is specified spline by the template argument <TT>ORDER</tt>.
        */
    explicit BSplineBase(unsigned int derivativeOrder = 0)
    : s1_(derivativeOrder)
    {}

        /** Unary function call.
            Returns the value of the spline with the derivative order given in the
            constructor. Note that only derivatives up to <tt>ORDER-1</tt> are
            continuous, and derivatives above <tt>ORDER+1</tt> are zero.
        */
    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder());
    }

        /** Binary function call.
            The given derivative order is added to the derivative order
            specified in the constructor. Note that only derivatives up to <tt>ORDER-1</tt> are
            continuous, and derivatives above <tt>ORDER+1</tt> are zero.
        */
    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder() + derivative_order);
    }

        /** Index operator. Same as unary function call.
        */
    value_type operator[](value_type x) const
        { return operator()(x); }

        /** Get the required filter radius for a discrete approximation of the
            spline. Always equal to <tt>(ORDER + 1) / 2.0</tt>.
        */
    double radius() const
        { return (ORDER + 1) * 0.5; }

        /** Get the derivative order of the Gaussian.
        */
    unsigned int derivativeOrder() const
        { return s1_.derivativeOrder(); }

        /** Get the prefilter coefficients required for interpolation.
            To interpolate with a B-spline, \ref resamplingConvolveImage()
            can be used. However, the image to be interpolated must be
            pre-filtered using \ref recursiveFilterX() and \ref recursiveFilterY()
            with the filter coefficients given by this function. The length of the array
            corresponds to how many times the above recursive filtering
            has to be applied (zero length means no prefiltering necessary).
        */
    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef ArrayVector<ArrayVector<T> > WeightMatrix;

        /** Get the coefficients to transform spline coefficients into
            the coefficients of the corresponding polynomial.
            Currently internally used in SplineImageView; needs more
            documentation ???
        */
    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;
    
        // factory function for the prefilter coefficients array
    static ArrayVector<double> calculatePrefilterCoefficients();

        // factory function for the weight matrix
    static WeightMatrix calculateWeightMatrix();

    BSplineBase<ORDER-1, T> s1_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <int ORDER, class T>
ArrayVector<double> BSplineBase<ORDER, T>::prefilterCoefficients_(BSplineBase<ORDER, T>::calculatePrefilterCoefficients());

template <int ORDER, class T>
typename BSplineBase<ORDER, T>::WeightMatrix BSplineBase<ORDER, T>::weightMatrix_(calculateWeightMatrix());

template <int ORDER, class T>
typename BSplineBase<ORDER, T>::result_type
BSplineBase<ORDER, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    if(derivative_order == 0)
    {
        T n12 = (ORDER + 1.0) / 2.0;
        return ((n12 + x) * s1_(x + 0.5) + (n12 - x) * s1_(x - 0.5)) / ORDER;
    }
    else
    {
        --derivative_order;
        return s1_(x + 0.5, derivative_order) - s1_(x - 0.5, derivative_order);
    }
}

template <int ORDER, class T>
ArrayVector<double> 
BSplineBase<ORDER, T>::calculatePrefilterCoefficients()
{
    ArrayVector<double> res;
    if(ORDER > 1)
    {
        const int r = ORDER / 2;
        StaticPolynomial<2*r, double> p(2*r);
        BSplineBase spline;
        for(int i = 0; i <= 2*r; ++i)
            p[i] = spline(T(i-r));
        ArrayVector<double> roots;
        polynomialRealRoots(p, roots);
        for(unsigned int i = 0; i < roots.size(); ++i)
            if(VIGRA_CSTD::fabs(roots[i]) < 1.0)
                res.push_back(roots[i]);
    }
    return res;
}

template <int ORDER, class T>
typename BSplineBase<ORDER, T>::WeightMatrix
BSplineBase<ORDER, T>::calculateWeightMatrix()
{
    WeightMatrix res(ORDER+1, ArrayVector<T>(ORDER+1));
    double faculty = 1.0;
    for(int d = 0; d <= ORDER; ++d)
    {
        if(d > 1)
            faculty *= d;
        double x = ORDER / 2; // (note: integer division)
        BSplineBase spline;
        for(int i = 0; i <= ORDER; ++i, --x)
            res[d][i] = spline(x, d) / faculty;
    }
    return res;
}

/********************************************************/
/*                                                      */
/*                     BSpline<N, T>                    */
/*                                                      */
/********************************************************/

/** Spline functors for arbitrary orders.

    Provides the interface of \ref vigra::BSplineBase with a more convenient
    name -- see there for more documentation.
*/
template <int ORDER, class T = double>
class BSpline
: public BSplineBase<ORDER, T>
{
  public:
        /** Constructor forwarded to the base class constructor..
        */
    explicit BSpline(unsigned int derivativeOrder = 0)
    : BSplineBase<ORDER, T>(derivativeOrder)
    {}
};

/********************************************************/
/*                                                      */
/*                     BSpline<0, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSplineBase<0, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 0 };

    explicit BSplineBase(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
         return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        return x.value < Value::ONE_HALF && -Value::ONE_HALF <= x.value
                   ? Value(Value::ONE, FPNoShift)
                   : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> const & x) const
    {
        return x < 0.5 && -0.5 <= x 
                   ? autodiff::DualVector<U, N>(1.0)
                   : autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 0.5; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[1][1];
  
    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const
    {
        if(derivative_order == 0)
            return x < 0.5 && -0.5 <= x ?
                     1.0
                   : 0.0;
        else
            return 0.0;
    }

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSplineBase<0, T>::prefilterCoefficients_;

template <class T>
typename BSplineBase<0, T>::WeightMatrix BSplineBase<0, T>::weightMatrix_ = {{ 1.0 }};

/********************************************************/
/*                                                      */
/*                     BSpline<1, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<1, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum  StaticOrder { order = 1 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        int v = abs(x.value);
        return v < Value::ONE ?
                Value(Value::ONE - v, FPNoShift)
                : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        return x < 1.0 
                    ? 1.0 - x
                    : autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 1.0; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[2][2];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    T exec(T x, unsigned int derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<1, T>::prefilterCoefficients_;

template <class T>
typename BSpline<1, T>::WeightMatrix BSpline<1, T>::weightMatrix_ = {{ 1.0, 0.0}, {-1.0, 1.0}};

template <class T>
T BSpline<1, T>::exec(T x, unsigned int derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                    1.0 - x
                    : 0.0;
        }
        case 1:
        {
            return x < 0.0 ?
                     -1.0 <= x ?
                          1.0
                     : 0.0
                   : x < 1.0 ?
                       -1.0
                     : 0.0;
        }
        default:
            return 0.0;
    }
}

/********************************************************/
/*                                                      */
/*                     BSpline<2, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<2, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 2 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        enum { ONE_HALF = Value::ONE_HALF, THREE_HALVES = ONE_HALF * 3, THREE_QUARTERS = THREE_HALVES / 2,
               PREMULTIPLY_SHIFT1 = FracBits <= 16 ? 0 : FracBits - 16,
               PREMULTIPLY_SHIFT2 = FracBits - 1 <= 16 ? 0 : FracBits - 17,
               POSTMULTIPLY_SHIFT1 = FracBits - 2*PREMULTIPLY_SHIFT1,
               POSTMULTIPLY_SHIFT2 = FracBits - 2*PREMULTIPLY_SHIFT2  };
        int v = abs(x.value);
        return v == ONE_HALF
                   ? Value(ONE_HALF, FPNoShift)
                   : v <= ONE_HALF
                       ? Value(THREE_QUARTERS -
                               (int)(sq((unsigned)v >> PREMULTIPLY_SHIFT2) >> POSTMULTIPLY_SHIFT2), FPNoShift)
                       : v < THREE_HALVES
                            ? Value((int)(sq((unsigned)(THREE_HALVES-v) >> PREMULTIPLY_SHIFT1) >> (POSTMULTIPLY_SHIFT1 + 1)), FPNoShift)
                            : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        return x < 0.5 
                   ? 0.75 - x*x
                   : x < 1.5 
                        ? 0.5 * sq(1.5 - x)
                        : autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 1.5; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[3][3];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<2, T>::prefilterCoefficients_(1, 2.0*M_SQRT2 - 3.0);

template <class T>
typename BSpline<2, T>::WeightMatrix BSpline<2, T>::weightMatrix_ = 
                           {{ 0.125, 0.75, 0.125},
                            {-0.5, 0.0, 0.5},
                            { 0.5, -1.0, 0.5}};

template <class T>
typename BSpline<2, T>::result_type
BSpline<2, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            return x < 0.5 ?
                    0.75 - x*x
                    : x < 1.5 ?
                        0.5 * sq(1.5 - x)
                    : 0.0;
        }
        case 1:
        {
            return x >= -0.5 ?
                     x <= 0.5 ?
                       -2.0 * x
                     : x < 1.5 ?
                         x - 1.5
                       : 0.0
                   : x > -1.5 ?
                       x + 1.5
                     : 0.0;
        }
        case 2:
        {
            return x >= -0.5 ?
                     x < 0.5 ?
                         -2.0
                     : x < 1.5 ?
                         1.0
                       : 0.0
                   : x >= -1.5 ?
                       1.0
                     : 0.0;
        }
        default:
            return 0.0;
    }
}

/********************************************************/
/*                                                      */
/*                     BSpline<3, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<3, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 3 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        enum { ONE = Value::ONE, TWO = 2 * ONE, TWO_THIRDS = TWO / 3, ONE_SIXTH = ONE / 6,
               PREMULTIPLY_SHIFT = FracBits <= 16 ? 0 : FracBits - 16,
               POSTMULTIPLY_SHIFT = FracBits - 2*PREMULTIPLY_SHIFT };
        int v = abs(x.value);
        return v == ONE
                   ? Value(ONE_SIXTH, FPNoShift)
                   : v < ONE
                       ? Value(TWO_THIRDS +
                               (((int)(sq((unsigned)v >> PREMULTIPLY_SHIFT) >> (POSTMULTIPLY_SHIFT + PREMULTIPLY_SHIFT))
                                       * (((v >> 1) - ONE) >> PREMULTIPLY_SHIFT)) >> POSTMULTIPLY_SHIFT), FPNoShift)
                       : v < TWO
                            ? Value((int)((sq((unsigned)(TWO-v) >> PREMULTIPLY_SHIFT) >> (POSTMULTIPLY_SHIFT + PREMULTIPLY_SHIFT))
                                      * ((unsigned)(TWO-v) >> PREMULTIPLY_SHIFT) / 6) >> POSTMULTIPLY_SHIFT, FPNoShift)
                            : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        if(x < 1.0)
        {
            return 2.0/3.0 + x*x*(-1.0 + 0.5*x);
        }
        else if(x < 2.0)
        {
            x = 2.0 - x;
            return x*x*x/6.0;
        }
        else
            return autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    result_type dxx(argument_type x) const
        { return operator()(x, 2); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 2.0; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[4][4];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<3, T>::prefilterCoefficients_(1, VIGRA_CSTD::sqrt(3.0) - 2.0);

template <class T>
typename BSpline<3, T>::WeightMatrix BSpline<3, T>::weightMatrix_ = 
                           {{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 0.0},
                            {-0.5, 0.0, 0.5, 0.0},
                            { 0.5, -1.0, 0.5, 0.0},
                            {-1.0 / 6.0, 0.5, -0.5, 1.0 / 6.0}};

template <class T>
typename BSpline<3, T>::result_type
BSpline<3, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x < 1.0)
            {
                return 2.0/3.0 + x*x*(-1.0 + 0.5*x);
            }
            else if(x < 2.0)
            {
                x = 2.0 - x;
                return x*x*x/6.0;
            }
            else
                return 0.0;
        }
        case 1:
        {
            double s = x < 0.0 ?
                         -1.0
                       :  1.0;
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                     s*x*(-2.0 + 1.5*x)
                   : x < 2.0 ?
                       -0.5*s*sq(2.0 - x)
                     : 0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                     3.0*x - 2.0
                   : x < 2.0 ?
                       2.0 - x
                     : 0.0;
        }
        case 3:
        {
            return x < 0.0 ?
                     x < -1.0 ?
                       x < -2.0 ?
                         0.0
                       : 1.0
                     : -3.0
                   : x < 1.0 ?
                       3.0
                     : x < 2.0 ?
                         -1.0
                       : 0.0;
        }
        default:
            return 0.0;
    }
}

typedef BSpline<3, double> CubicBSplineKernel;

/********************************************************/
/*                                                      */
/*                     BSpline<4, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<4, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 4 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        if(x <= 0.5)
        {
            return 115.0/192.0 + x*x*(-0.625 + x*x*0.25);
        }
        else if(x < 1.5)
        {
            return (55.0/16.0 + x*(1.25 + x*(-7.5 + x*(5.0 - x)))) / 6.0;
        }
        else if(x < 2.5)
        {
            x = 2.5 - x;
            return sq(x*x) / 24.0;
        }
        else
            return autodiff::DualVector<U, N>(0.0);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    result_type dxx(argument_type x) const
        { return operator()(x, 2); }

    result_type dx3(argument_type x) const
        { return operator()(x, 3); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 2.5; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[5][5];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    static ArrayVector<double> calculatePrefilterCoefficients()
    {
        ArrayVector<double> b(2);
        // -19 + 4*sqrt(19) + 2*sqrt(2*(83 - 19*sqrt(19)))
        b[0] = -0.361341225900220177092212841325;
        // -19 - 4*sqrt(19) + 2*sqrt(2*(83 + 19*sqrt(19)))
        b[1] = -0.013725429297339121360331226939;
        return b;
    }

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<4, T>::prefilterCoefficients_(calculatePrefilterCoefficients());

template <class T>
typename BSpline<4, T>::WeightMatrix BSpline<4, T>::weightMatrix_ = 
                           {{ 1.0/384.0, 19.0/96.0, 115.0/192.0, 19.0/96.0, 1.0/384.0},
                            {-1.0/48.0, -11.0/24.0, 0.0, 11.0/24.0, 1.0/48.0},
                            { 1.0/16.0, 1.0/4.0, -5.0/8.0, 1.0/4.0, 1.0/16.0},
                            {-1.0/12.0, 1.0/6.0, 0.0, -1.0/6.0, 1.0/12.0},
                            { 1.0/24.0, -1.0/6.0, 0.25, -1.0/6.0, 1.0/24.0}};

template <class T>
typename BSpline<4, T>::result_type
BSpline<4, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return 115.0/192.0 + x*x*(-0.625 + x*x*0.25);
            }
            else if(x < 1.5)
            {
                return (55.0/16.0 + x*(1.25 + x*(-7.5 + x*(5.0 - x)))) / 6.0;
            }
            else if(x < 2.5)
            {
                x = 2.5 - x;
                return sq(x*x) / 24.0;
            }
            else
                return 0.0;
        }
        case 1:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return s*x*(-1.25 + x*x);
            }
            else if(x < 1.5)
            {
                return s*(5.0 + x*(-60.0 + x*(60.0 - 16.0*x))) / 24.0;
            }
            else if(x < 2.5)
            {
                x = 2.5 - x;
                return s*x*x*x / -6.0;
            }
            else
                return 0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return -1.25 + 3.0*x*x;
            }
            else if(x < 1.5)
            {
                return -2.5 + x*(5.0 - 2.0*x);
            }
            else if(x < 2.5)
            {
                x = 2.5 - x;
                return x*x / 2.0;
            }
            else
                return 0.0;
        }
        case 3:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return s*x*6.0;
            }
            else if(x < 1.5)
            {
                return s*(5.0 - 4.0*x);
            }
            else if(x < 2.5)
            {
                return s*(x - 2.5);
            }
            else
                return 0.0;
        }
        case 4:
        {
            return x < 0.0
                     ? x < -2.5
                         ? 0.0
                         : x < -1.5
                             ? 1.0
                             : x < -0.5
                                 ? -4.0
                                 : 6.0
                     : x < 0.5 
                         ? 6.0
                         : x < 1.5
                             ? -4.0
                             : x < 2.5
                                 ? 1.0
                                 : 0.0;
        }
        default:
            return 0.0;
    }
}

typedef BSpline<4, double> QuarticBSplineKernel;

/********************************************************/
/*                                                      */
/*                     BSpline<5, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<5, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 5 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        if(x <= 1.0)
        {
            return 0.55 + x*x*(-0.5 + x*x*(0.25 - x/12.0));
        }
        else if(x < 2.0)
        {
            return 17.0/40.0 + x*(0.625 + x*(-1.75 + x*(1.25 + x*(-0.375 + x/24.0))));
        }
        else if(x < 3.0)
        {
            x = 3.0 - x;
            return x*sq(x*x) / 120.0;
        }
        else
            return autodiff::DualVector<U, N>(0.0);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    result_type dxx(argument_type x) const
        { return operator()(x, 2); }

    result_type dx3(argument_type x) const
        { return operator()(x, 3); }

    result_type dx4(argument_type x) const
        { return operator()(x, 4); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 3.0; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[6][6];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    static ArrayVector<double> calculatePrefilterCoefficients()
    {
        ArrayVector<double> b(2);
        // -(13/2) + sqrt(105)/2 + sqrt(1/2*((135 - 13*sqrt(105))))
        b[0] = -0.430575347099973791851434783493;
        // (1/2)*((-13) - sqrt(105) + sqrt(2*((135 + 13*sqrt(105)))))
        b[1] = -0.043096288203264653822712376822;
        return b;
    }

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<5, T>::prefilterCoefficients_(calculatePrefilterCoefficients());

template <class T>
typename BSpline<5, T>::WeightMatrix BSpline<5, T>::weightMatrix_ = 
                           {{ 1.0/120.0, 13.0/60.0, 11.0/20.0, 13.0/60.0, 1.0/120.0, 0.0},
                            {-1.0/24.0, -5.0/12.0, 0.0, 5.0/12.0, 1.0/24.0, 0.0},
                            { 1.0/12.0, 1.0/6.0, -0.5, 1.0/6.0, 1.0/12.0, 0.0},
                            {-1.0/12.0, 1.0/6.0, 0.0, -1.0/6.0, 1.0/12.0, 0.0},
                            { 1.0/24.0, -1.0/6.0, 0.25, -1.0/6.0, 1.0/24.0, 0.0},
                            {-1.0/120.0, 1.0/24.0, -1.0/12.0, 1.0/12.0, -1.0/24.0, 1.0/120.0}};

template <class T>
typename BSpline<5, T>::result_type
BSpline<5, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return 0.55 + x*x*(-0.5 + x*x*(0.25 - x/12.0));
            }
            else if(x < 2.0)
            {
                return 17.0/40.0 + x*(0.625 + x*(-1.75 + x*(1.25 + x*(-0.375 + x/24.0))));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return x*sq(x*x) / 120.0;
            }
            else
                return 0.0;
        }
        case 1:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return s*x*(-1.0 + x*x*(1.0 - 5.0/12.0*x));
            }
            else if(x < 2.0)
            {
                return s*(0.625 + x*(-3.5 + x*(3.75 + x*(-1.5 + 5.0/24.0*x))));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return s*sq(x*x) / -24.0;
            }
            else
                return 0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return -1.0 + x*x*(3.0 -5.0/3.0*x);
            }
            else if(x < 2.0)
            {
                return -3.5 + x*(7.5 + x*(-4.5 + 5.0/6.0*x));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return x*x*x / 6.0;
            }
            else
                return 0.0;
        }
        case 3:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return s*x*(6.0 - 5.0*x);
            }
            else if(x < 2.0)
            {
                return s*(7.5 + x*(-9.0 + 2.5*x));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return -0.5*s*x*x;
            }
            else
                return 0.0;
        }
        case 4:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return 6.0 - 10.0*x;
            }
            else if(x < 2.0)
            {
                return -9.0 + 5.0*x;
            }
            else if(x < 3.0)
            {
                return 3.0 - x;
            }
            else
                return 0.0;
        }
        case 5:
        {
            return x < 0.0 ?
                     x < -2.0 ?
                       x < -3.0 ?
                         0.0
                       : 1.0
                     : x < -1.0 ?
                         -5.0
                       : 10.0
                   : x < 2.0 ?
                       x < 1.0 ?
                         -10.0
                       : 5.0
                     : x < 3.0 ?
                         -1.0
                       : 0.0;
        }
        default:
            return 0.0;
    }
}

typedef BSpline<5, double> QuinticBSplineKernel;

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

/********************************************************/
/*                                                      */
/*                      CatmullRomSpline                */
/*                                                      */
/********************************************************/

/** Interpolating 3-rd order splines.

    Implements the Catmull/Rom cardinal function

    \f[ f(x) = \left\{ \begin{array}{ll}
                                  \frac{3}{2}x^3 - \frac{5}{2}x^2 + 1 & |x| \leq 1 \\
                                  -\frac{1}{2}x^3 + \frac{5}{2}x^2 -4x + 2 & |x| \leq 2 \\
                                  0 & \mbox{otherwise}
                        \end{array}\right.
    \f]

    It can be used as a functor, and as a kernel for
    \ref resamplingConvolveImage() to create a differentiable interpolant
    of an image. However, it should be noted that a twice differentiable
    interpolant can be created with only slightly more effort by recursive
    prefiltering followed by convolution with a 3rd order B-spline.

    <b>\#include</b> \<vigra/splines.hxx\><br>
    Namespace: vigra
*/
template <class T = double>
class CatmullRomSpline
{
public:
        /** the kernel's value type
        */
    typedef T value_type;
        /** the unary functor's argument type
        */
    typedef T argument_type;
        /** the unary functor's result type
        */
    typedef T result_type;
        /** the splines polynomial order
        */
    enum StaticOrder { order = 3 };

        /** function (functor) call
        */
    result_type operator()(argument_type x) const;

        /** index operator -- same as operator()
        */
    T operator[] (T x) const
        { return operator()(x); }

        /** Radius of the function's support.
            Needed for  \ref resamplingConvolveImage(), always 2.
        */
    int radius() const
        {return 2;}

        /** Derivative order of the function: always 0.
        */
    unsigned int derivativeOrder() const
        { return 0; }

        /** Prefilter coefficients for compatibility with \ref vigra::BSpline.
            (array has zero length, since prefiltering is not necessary).
        */
    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }
    
  protected:
    static ArrayVector<double> prefilterCoefficients_;
};

template <class T>
ArrayVector<double> CatmullRomSpline<T>::prefilterCoefficients_;

template <class T>
typename CatmullRomSpline<T>::result_type
CatmullRomSpline<T>::operator()(argument_type x) const
{
    x = VIGRA_CSTD::fabs(x);
    if (x <= 1.0)
    {
        return 1.0 + x * x * (-2.5 + 1.5 * x);
    }
    else if (x >= 2.0)
    {
        return 0.0;
    }
    else
    {
        return 2.0 + x * (-4.0 + x * (2.5 - 0.5 * x));
    }
}


//@}


//@{

	namespace RBF {

		namespace Kernels {

		    template<int DIM> 
		    struct TPSKernel {
		    };
		    template<>
		    struct TPSKernel<1> {
			template<class VALUETYPE>
			static VALUETYPE eval(const VALUETYPE &r) {
			    return r * r * r;
			}
		    };
#if 0
		    template<>
		    struct TPSKernel<2> {
			template<class VALUETYPE>
			static VALUETYPE operator()(const VALUETYPE &r) {
				if(r > 0.0)
					return r * r * log(r);
				else
					return 0.0;
			}
		    };
		    template<>
		    struct TPSKernel<3> {
			template<class VALUETYPE>
			static VALUETYPE operator()(const VALUETYPE &r) {
			    return r;
			}
		    };
#endif
		} // namespace Kernels


		//! Simple class to compute a thin-plate spline (TPS) mapping
		// (approximation or interpolation, depending on @param relaxation)
		template<int INPUTDIM, int OUTPUTDIM, class VALUETYPE>
		class ThinPlateSpline {
		public:
			typedef MultiArray<2, VALUETYPE> Matrix;
			typedef MultiArray<2, VALUETYPE> KnotMatrix;
			typedef MultiArray<2, VALUETYPE> ValueMatrix;
			typedef MultiArrayView<2, VALUETYPE> KnotMatrixView;
			typedef MultiArrayView<2, VALUETYPE> ValueMatrixView;
 		        typedef MultiArray<1, VALUETYPE> Vector;
 		        typedef MultiArrayView<1, VALUETYPE, StridedArrayTag> VectorView;
			typedef typename MultiArray<1, VALUETYPE>::difference_type C1;
			typedef typename MultiArray<2, VALUETYPE>::difference_type C2;

		    typedef Kernels::TPSKernel<INPUTDIM> Kernel;

		    ThinPlateSpline(const KnotMatrixView &knots, const ValueMatrixView &targetValues, const VALUETYPE &relaxation = 0.0) : m_Knots(knots) 
			{
			    vigra_precondition(relaxation >= 0.0, "relaxation parameter must be nonnegative");

				vigra_precondition(INPUTDIM >= 1, "unsupported input dimension");
				vigra_precondition(INPUTDIM <= 3, "unsupported input dimension");
				vigra_precondition(OUTPUTDIM >= 1, "unsupported output dimension");
				vigra_precondition(OUTPUTDIM <= 3, "unsupported output dimension");
				vigra_precondition(knots.shape(1) == INPUTDIM, "incompatible knot array dimension");
				vigra_precondition(targetValues.shape(1) == OUTPUTDIM, "incompatible target value dimension");

				const int Npoints = knots.shape(0);
				std::cout << " knots shape: " << knots.shape(0) << "," << knots.shape(1) << std::endl;
				std::cout << " targetValues shape: " << targetValues.shape(0) << "," << targetValues.shape(1) << std::endl;
				vigra_precondition(targetValues.shape(0) == Npoints, "incompatible knot and target arrays");
				const int affineOffset = 1 + INPUTDIM;
				const int tpsdim = affineOffset + Npoints;
				
				Matrix systemMatrix(C2(tpsdim, tpsdim), 0.0);
				// Setup the linear system:
				// (0  P)   (cp)   (0)
				// (Pt K) . (ck) = (b)

				// affine components go first: constant, x1, ...; then follow the kernel pts.
				// fill system matrix by columns
				for (int i=0; i < Npoints; ++i) {
				    EvalBasisForPoint(knots.template bind<0>(i), systemMatrix.template bind<1>(affineOffset + i));
				}
				// copy the polynomial part of the system matrix (transpose)
				systemMatrix.subarray(C2(affineOffset, 0), C2(tpsdim, affineOffset)).copy(
				    systemMatrix.subarray(C2(0, affineOffset), C2(affineOffset, tpsdim)).transpose()
				    );
				

				Matrix rhs(C2(tpsdim, OUTPUTDIM), 0.0);
				rhs.subarray(C2(affineOffset, 0), C2(tpsdim, OUTPUTDIM)).copy(targetValues);

				{
				    for (int i=0; i < tpsdim; ++i) {
					for (int j=0; j < tpsdim; ++j) {
					    std::cout << "\t" << systemMatrix(i, j);
					}
					std::cout << "   | "; 
					for (int j=0; j < OUTPUTDIM; ++j) {
					    std::cout << "\t" << rhs(i, j);
					}
					std::cout << std::endl;
				    }
				}


				m_Coeffs.reshape(C2(tpsdim, OUTPUTDIM), 0.0);

				// compute QR decomposition
				const int rank = linalg::linearSolveQRReplace(systemMatrix, rhs, m_Coeffs);
				if (rank != tpsdim) {
				    std::cout << " QR solve return rank=" << rank << std::endl;
				}
				std::cout << "coeffs: ";
				std::copy(m_Coeffs.begin(), m_Coeffs.end(), std::ostream_iterator<double>(std::cout, ","));
				std::cout << std::endl;
			}

		        void EvalBasisForPoint(const VectorView& point, VectorView result) const {
			    const int Npoints = m_Knots.shape(0);
			    const int affineOffset = 1 + INPUTDIM;
			    const int tpsdim = affineOffset + Npoints;
			    
			    vigra_precondition(point.shape(0) == INPUTDIM, "wrong input dimension");
			    vigra_precondition(result.shape(0) == tpsdim, "wrong output dimension");
			    
			    result(0) = 1.0;
			    for (int i = 0; i < INPUTDIM; ++i) {
				result(1+i) = point(i);
			    }
			    // kernel evaluation
			    for (int i = 0; i < Npoints; ++i) {
				const VectorView kernelPoint = m_Knots.template bind<0>(i);
				Vector diff(kernelPoint);
				diff -= point;
				const double norm = sqrt(dot(diff, diff));
				result(affineOffset + i) = Kernel::eval(norm >= 0 ? norm : 0.0);
			    }			    
			}
		    
		    const int tpsDim() const {
			    const int Npoints = m_Knots.shape(0);
			    const int affineOffset = 1 + INPUTDIM;
			    return affineOffset + Npoints;
		    }

			bool operator()(const KnotMatrixView &evalPositions, ValueMatrixView &resultValues) const {
				const int Npoints = evalPositions.shape(0);
				vigra_precondition(evalPositions.shape(1) == INPUTDIM, "eval position dimension mismatch");
				vigra_precondition(resultValues.shape(0) == Npoints, "result array shape does not match");
				vigra_precondition(resultValues.shape(1) == OUTPUTDIM, "result array shape does not match");

				const int tpsdim = tpsDim();

				Vector tempVector(C1(tpsdim), 0.0);

				for (int i=0; i < Npoints; ++i) {
				    EvalBasisForPoint(evalPositions.template bind<0>(i), tempVector);
				    for (int j=0; j < OUTPUTDIM; ++j) {
					resultValues(i,j) = dot(m_Coeffs.template bind<1>(j), tempVector);
				    }
				}
				return true;
			}

			ValueMatrix operator()(const KnotMatrixView &evalPositions) const {
				const int Npoints = evalPositions.shape(0);
				ValueMatrix result(C2(Npoints, OUTPUTDIM), 0.0);
				if (this->operator()(evalPositions, result))
					return result;
				else 
				    throw std::runtime_error("TPS evaluation failed.");
			}

		protected:
			void EvaluateTPSBasis() {}
			void ComputeInterpolationMatrix() {}

		protected:
			Matrix m_Coeffs;
			KnotMatrix m_Knots;
		};


	} // namespace RBF

/*
	template<class T1, class T2, int DIM>
	double
	RBF(const MultiArrayView<2, T1> & v1, const MultiArrayView<2, T2> & v2)
	{
		boost::numeric::ublas::vector<double> v(Dim);
		for(unsigned int i(0); i < v.size(); ++i)
		{
			v(i) = v1[i] - v2[i];
		}
		return radial_basis(norm_2(v));
	}
*/


//@}



} // namespace vigra


#endif /* VIGRA_SPLINES_HXX */
