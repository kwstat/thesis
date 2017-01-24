// IntervalComplex.h
//   Kevin Wright

#ifndef __INTCOMPLEX__
#define __INTCOMPLEX__

#include "Functions.h"  // For Sqr(INTERVAL)
#include "Interval.h"
#include "Complex.h"
#include "Portab.h"
#include "BiasF.h"  // For RoundUp, RoundDown
#include "Configuration.h"
#include <math.h>
#include <iostream.h>

class INTERVAL_COMPLEX {
  INTERVAL re;
  INTERVAL im;
public:
  // Constructor
  INTERVAL_COMPLEX (){}

  // Various copy constructors
  INTERVAL_COMPLEX (REAL r){
    re = Hull(r); 
    im = Hull(0.0);
  }
  INTERVAL_COMPLEX (REAL r, REAL i){
    re = Hull(r); 
    im = Hull(i);
  }
  INTERVAL_COMPLEX (INTERVAL r){
    re = r; 
    im = Hull(0.00);
  }
  INTERVAL_COMPLEX (INTERVAL r, INTERVAL i){
    re = r; 
    im = i;
  }
  INTERVAL_COMPLEX (COMPLEX z){
    re = Hull(Re(z));
    im = Hull(Im(z));
	}

  // No deconstructor is specified

  // First functions for REAL types
  INTERVAL_COMPLEX & operator= (const REAL& x){
    re = Hull(x); 
    im = Hull(0.0); 
    return *this; }
  INTERVAL_COMPLEX & operator+= (const REAL& x){
    re += x;
    return *this; }
  INTERVAL_COMPLEX & operator-= (const REAL& x){
    re -= x;
    return *this; }
  INTERVAL_COMPLEX & operator*= (const REAL& x){
		re *= x;
    im *= x;
    return *this; }
  INTERVAL_COMPLEX & operator/= (const REAL& x){
    re /= x;
    im /= x;
    return *this; }
    
  // Then functions for INTERVAL class
  INTERVAL_COMPLEX & operator= (const INTERVAL& x){
    re = x;
    im = Hull(0.0);
    return *this; }
  INTERVAL_COMPLEX & operator+= (const INTERVAL& x){
    re += x;
    return *this;}
  INTERVAL_COMPLEX & operator-= (const INTERVAL& x){
    re -= x;
    return *this;}
  INTERVAL_COMPLEX & operator*= (const INTERVAL& x){
    re *= x;
    im *= x; 
    return *this;}
  INTERVAL_COMPLEX & operator/= (const INTERVAL& x){
    re /= x;
    im /= x;
    return *this; }
// Note, the above definition assumes 0 not in X


  // Functions for COMPLEX type
  INTERVAL_COMPLEX & operator= (COMPLEX &x)
		{ re = Re(x); im = Im(x); return *this; }
  INTERVAL_COMPLEX & operator+= (COMPLEX &x)
    { re += Re(x); im += Im(x); return *this; }
  INTERVAL_COMPLEX & operator-= (COMPLEX &x)
    { re -= Re(x); im -= Im(x); return *this; }
  INTERVAL_COMPLEX & operator*= (COMPLEX &x)
    { INTERVAL t = re;
      re = re * Re(x) - im * Im(x);
      im = im * Re(x) + t  * Im(x);
      return *this;
    }
  INTERVAL_COMPLEX & operator/= (COMPLEX &x);


  // Finally functions for INTERVAL_COMPLEX
  INTERVAL_COMPLEX & operator= (const INTERVAL_COMPLEX& x)
		{ re = Re(x); im = Im(x); return *this; }
  INTERVAL_COMPLEX & operator+= (const INTERVAL_COMPLEX& x)
    { re += Re(x); im += Im(x); return *this; }
  INTERVAL_COMPLEX & operator-= (const INTERVAL_COMPLEX& x)
    { re -= Re(x); im -= Im(x); return *this; }
  INTERVAL_COMPLEX & operator*= (const INTERVAL_COMPLEX& x)
    { INTERVAL t = re;
      re = re * Re(x) - im * Im(x);
      im = im * Re(x) + t  * Im(x);
      return *this;
    }
  INTERVAL_COMPLEX & operator/= (INTERVAL_COMPLEX& x);

  // Unary operators
  friend INTERVAL_COMPLEX operator+ (const INTERVAL_COMPLEX& x) { return x; }
  friend INTERVAL_COMPLEX operator- (const INTERVAL_COMPLEX& x)
    { return INTERVAL_COMPLEX (-Re(x), -Im(x)); }

  friend INTERVAL Re (const INTERVAL_COMPLEX &x) { return x.re; }
  friend INTERVAL Im (const INTERVAL_COMPLEX &x) { return x.im; }

  // Binary operators
  // Addition----------------------------------------------
  // INTERVAL_COMPLEX & REAL
  friend INTERVAL_COMPLEX operator+ (const INTERVAL_COMPLEX& x, REAL y){
    return INTERVAL_COMPLEX (Re(x) + y, Im(x)); }
  friend INTERVAL_COMPLEX operator+ (REAL x, const INTERVAL_COMPLEX& y){
     return INTERVAL_COMPLEX (x + Re(y), Im(y)); }
  // INTERVAL_COMPLEX & INTERVAL
  friend INTERVAL_COMPLEX operator+ (const INTERVAL_COMPLEX& x, INTERVAL y){
    return INTERVAL_COMPLEX (Re(x) + y, Im(x)); }
  friend INTERVAL_COMPLEX operator+ (INTERVAL x, const INTERVAL_COMPLEX& y){
     return INTERVAL_COMPLEX (x + Re(y), Im(y)); }  
  // INTERVAL_COMPLEX  &  COMPLEX 
  friend INTERVAL_COMPLEX operator+ (const INTERVAL_COMPLEX& x, 
																		 COMPLEX& y){ 
    return INTERVAL_COMPLEX (Re(x) + Re(y), Im(x)+Im(y)); }
  friend INTERVAL_COMPLEX operator+ (COMPLEX& x, 
																		 const INTERVAL_COMPLEX& y){ 
		return INTERVAL_COMPLEX (Re(x) + Re(y), Im(x) + Im(y)); }
  // INTERVAL_COMPLEX & INTERVAL_COMPLEX
  friend INTERVAL_COMPLEX operator+ (const INTERVAL_COMPLEX& x, 
																			const INTERVAL_COMPLEX& y){ 
		return INTERVAL_COMPLEX (Re(x) + Re(y), Im(x) + Im(y)); }

  // Subtraction-----------------------------------------------
  // INTERVAL_COMPLEX & REAL 
  friend INTERVAL_COMPLEX operator- (const INTERVAL_COMPLEX &x, REAL y)
    { return INTERVAL_COMPLEX (Re(x) - y, Im(x)); }
  friend INTERVAL_COMPLEX operator- (REAL x, const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (x - Re(y), -Im(y)); }
  // INTERVAL_COMPLEX & INTERVAL 
  friend INTERVAL_COMPLEX operator- (const INTERVAL_COMPLEX &x, INTERVAL y)
    { return INTERVAL_COMPLEX (Re(x) - y, Im(x)); }
  friend INTERVAL_COMPLEX operator- (INTERVAL x, const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (x - Re(y), -Im(y)); }
  // INTERVAL_COMPLEX & COMPLEX
  friend INTERVAL_COMPLEX operator- (const INTERVAL_COMPLEX& x, 
																		 COMPLEX& y)
    { return INTERVAL_COMPLEX (Re(x) - Re(y), Im(x)-Im(y)); }
  friend INTERVAL_COMPLEX operator- (COMPLEX& x, 
																		 const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (Re(x) - Re(y), Im(x)-Im(y)); }
  // INTERVAL_COMPLEX & INTERVAL_COMPLEX
  friend INTERVAL_COMPLEX operator- (const INTERVAL_COMPLEX& x, 
																		 const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (Re(x) - Re(y), Im(x) - Im(y)); }

  //  Multiplication ------------------------------------------------
  // INTERVAL_COMPLEX & REAL
  friend INTERVAL_COMPLEX operator* (const INTERVAL_COMPLEX& x, REAL y)
    { return INTERVAL_COMPLEX (Re(x) * y, Im(x) * y); }
  friend INTERVAL_COMPLEX operator* (REAL x, const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (x * Re(y), x * Im(y)); }
  // INTERVAL_COMPLEX & INTERVAL
  friend INTERVAL_COMPLEX operator* (const INTERVAL_COMPLEX& x, INTERVAL y)
    { return INTERVAL_COMPLEX (Re(x) * y, Im(x) * y); }
  friend INTERVAL_COMPLEX operator* (INTERVAL x, const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (x * Re(y), x * Im(y)); }
  // INTERVAL_COMPLEX & COMPLEX
  friend INTERVAL_COMPLEX operator* (const INTERVAL_COMPLEX& x, 
																		 COMPLEX& y)
    { return INTERVAL_COMPLEX (Re(x) * Re(y) - Im(x) * Im(y),
		      Im(x) * Re(y) + Im(y) * Re(x)); }
  friend INTERVAL_COMPLEX operator* (COMPLEX& x, 
																		 const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (Re(x) * Re(y) - Im(x) * Im(y),
		      Im(x) * Re(y) + Im(y) * Re(x)); }
  // INTERVAL_COMPLEX & INTERVAL_COMPLEX
  friend INTERVAL_COMPLEX operator* (const INTERVAL_COMPLEX& x, 
																		 const INTERVAL_COMPLEX& y)
    { return INTERVAL_COMPLEX (Re(x) * Re(y) - Im(x) * Im(y),
		      Im(x) * Re(y) + Im(y) * Re(x)); }

  // Division -----------------------------------------------------
  // INTERVAL_COMPLEX & REAL 
  friend INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, REAL y)
    { return INTERVAL_COMPLEX (Re(x) / y, Im(x) / y); }
  friend INTERVAL_COMPLEX operator/ (REAL x, const INTERVAL_COMPLEX& y);
  // INTERVAL_COMPLEX & INTERVAL
  friend INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, INTERVAL y)
    { return INTERVAL_COMPLEX (Re(x) / y, Im(x) / y); }
  friend INTERVAL_COMPLEX operator/ (INTERVAL x, const INTERVAL_COMPLEX& y);
  // INTERVAL_COMPLEX & COMPLEX
  friend INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, 
																		 COMPLEX& y);
  friend INTERVAL_COMPLEX operator/ (COMPLEX& x, 
																		 const INTERVAL_COMPLEX& y);
  // INTERVAL_COMPLEX & INTERVAL_COMPLEX 
  friend INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, 
																		 const INTERVAL_COMPLEX& y);

  // Logical comparisons
  // There's potential (but no need?) to expand these for different types
  friend INT operator== (const INTERVAL_COMPLEX& x, 
													const INTERVAL_COMPLEX& y)
    { return (Re(x) == Re(y)) && (Im(x) == Im(y)); }
  friend INT operator!= (const INTERVAL_COMPLEX& x, 
													const INTERVAL_COMPLEX& y)
    { return (Re(x) != Re(y)) || (Im(x) != Im(y)); }


  // Miscellany
  friend INTERVAL_COMPLEX Sqr (const INTERVAL_COMPLEX& x)
    { return INTERVAL_COMPLEX (Sqr(Re(x))-Sqr(Im(x)),2.0*Re(x)*Im(x)); }
  friend INTERVAL_COMPLEX Conjg (const INTERVAL_COMPLEX& x)
    { return INTERVAL_COMPLEX (Re(x), -Im(x)); }
  friend INTERVAL Abs (const INTERVAL_COMPLEX& x);


  // Stream operators.  Defined in *.C file
  friend ostream & operator<< (ostream &o, const INTERVAL_COMPLEX &x);
  //friend istream & operator>> (istream &i, INTERVAL_COMPLEX &x);

};  // end class

extern COMPLEX I;
#endif  /* __INTCOMPLEX__ */
