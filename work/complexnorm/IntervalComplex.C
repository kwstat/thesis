// IntervalComplex.C
// Kevin Wright

#include "IntervalComplex.h"
#include "Functions.h"  // For Sqr
#define NAIVE
//#define ROKNE
//#define DEBUG

// Denominator is COMPLEX.  Use the same ideas as in the PROFIL package.

// The following code gives an error:
// IntervalComplex.C:12: aggregate value used where a float was expected
// Line 12: REAL thisRe = re;
/*
INTERVAL_COMPLEX & INTERVAL_COMPLEX::operator /= (COMPLEX & x)
{
  REAL r, den;
  REAL thisRe = re;

  if (fabs (Re(x)) >= fabs (Im(x))) {
    r = Im(x) / Re(x);
    den = Re(x) + r * Im(x);
    re = (thisRe + r * im) / den;
    im = (im - r * thisRe) / den;
  }
  else {
    r = Re(x) / Im(x);
    den = Im(x) + r * Re(x);
    re = (r * thisRe + im) / den;
    im = (r * im - thisRe) / den;
  }
  return *this;
}
*/

INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, COMPLEX& y)
// PROFIL-style computation
{
  INTERVAL r,s,den;

  if (0 <= Re(y)  && !(0 <= Im(y))) {
    r = Re(y) /Im(y);
    den = Im(y) + r * Re(y);
    return INTERVAL_COMPLEX( (Re(x)*r + Im(x))/den, (Im(x)*r - Re(x))/den ); 
	}
  else if (0 <= Im(y) && !(0 <= Re(y))){
    r = Im(y)/Re(y);
    den = Re(y) + r * Im(y);
    return INTERVAL_COMPLEX( (Re(x) + r*Im(x))/den, (Im(x) - r*Re(x))/den );
	}
  else {
    r = Re(y) /Im(y);
		s = Im(y) /Re(y);
    if ( Sup(Abs(r)) < Sup(Abs(s)) ){
      den = Im(y) + r * Re(y);
      return INTERVAL_COMPLEX( (Re(x)*r + Im(x))/den, (Im(x)*r - Re(x))/den ); 
		}
		else {
      den = Re(y) + r * Im(y);
      return INTERVAL_COMPLEX( (Re(x) + r*Im(x))/den, (Im(x) - r*Re(x))/den );
		}
	}

}

#ifdef PROFIL
// Denominator is INTERVAL_COMPLEX
INTERVAL_COMPLEX & INTERVAL_COMPLEX::operator /= (INTERVAL_COMPLEX & x)
{
  INTERVAL r, den;
  INTERVAL thisRe = re;

  if (fabs (Re(x)) >= fabs (Im(x))) {
    r = Im(x) / Re(x);
    den = Re(x) + r * Im(x);
    re = (thisRe + r * im) / den;
    im = (im - r * thisRe) / den;
  }
  else {
    r = Re(x) / Im(x);
    den = Im(x) + r * Re(x);
    re = (r * thisRe + im) / den;
    im = (r * im - thisRe) / den;
  }
  return *this;
}
INTERVAL_COMPLEX operator/ (REAL x, const INTERVAL_COMPLEX& y)
{ 
  // THIS COULD BE IMPROVED
  /* fixme */
  INTERVAL den = Sqr(Re(y)) + Sqr(Im(y));
  return INTERVAL_COMPLEX (x*Re(y)/den, -x*Im(y)/den);
}
INTERVAL_COMPLEX operator/ (INTERVAL x, const INTERVAL_COMPLEX& y)
{
  // THIS COULD BE IMPROVED 
  /* fixme */
  INTERVAL den = Sqr(Re(y)) + Sqr(Im(y));
  return INTERVAL_COMPLEX (x*Re(y)/den, -x*Im(y)/den);
}
INTERVAL_COMPLEX operator/ (COMPLEX& x, const INTERVAL_COMPLEX& y)
{
  INTERVAL r,s,den;

  if (0 <= Re(y)  && !(0 <= Im(y))) {
    r = Re(y) /Im(y);
    den = Im(y) + r * Re(y);
    return INTERVAL_COMPLEX( (Re(x)*r + Im(x))/den, (Im(x)*r - Re(x))/den ); 
	}
  else if (0 <= Im(y) && !(0 <= Re(y))){
    r = Im(y)/Re(y);
    den = Re(y) + r * Im(y);
    return INTERVAL_COMPLEX( (Re(x) + r*Im(x))/den, (Im(x) - r*Re(x))/den );
	}
  else {
    r = Re(y) /Im(y);
		s = Im(y) /Re(y);
    if ( Sup(Abs(r)) < Sup(Abs(s)) ){
      den = Im(y) + r * Re(y);
      return INTERVAL_COMPLEX( (Re(x)*r + Im(x))/den, (Im(x)*r - Re(x))/den ); 
		}
		else {
      den = Re(y) + r * Im(y);
      return INTERVAL_COMPLEX( (Re(x) + r*Im(x))/den, (Im(x) - r*Re(x))/den );
		}
	}

}
#endif // PROFIL
//----------------------------------------------------------

#ifdef NAIVE
INTERVAL_COMPLEX & INTERVAL_COMPLEX::operator /= (INTERVAL_COMPLEX & x){
// /* fixme */
}
INTERVAL_COMPLEX operator/ (REAL x, const INTERVAL_COMPLEX& y){
  INTERVAL den = Sqr(Re(y)) + Sqr(Im(y));
  return INTERVAL_COMPLEX (x*Re(y)/den, -x*Im(y)/den);
}
INTERVAL_COMPLEX operator/ (INTERVAL x, const INTERVAL_COMPLEX& y){
  INTERVAL den = Sqr(Re(y)) + Sqr(Im(y));
  return INTERVAL_COMPLEX (x*Re(y)/den, -x*Im(y)/den);
}
INTERVAL_COMPLEX operator/ (COMPLEX& x, const INTERVAL_COMPLEX& y){
  INTERVAL den = Sqr(Re(y)) + Sqr(Im(y));
  return INTERVAL_COMPLEX ( (Re(x)*Re(y) + Im(x)*Im(y))/den,
													  (Im(x)*Re(y) - Re(x)*Im(y))/den );
}
INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, 
														const INTERVAL_COMPLEX& y){
  if (Abs(Inf(Re(x))) > 1.0e152){
    INTERVAL den = Sqr(Re(y)) + Sqr(Im(y));
/*
    cout << endl;
    cout << "Inside division: " << endl;
    cout << "x: " << x << endl;
    cout << "y: " << y << endl;
    cout << "Re(y): " << Re(y) << endl;
    cout << "Sqr(Re(y)): " << Sqr(Re(y)) << endl;
    cout << "Im(y): " << Im(y) << endl;
    cout << "Sqr(Im(y)): " << Sqr(Im(y)) << endl;
    cout << "Re(x)*Re(y): " << Re(x)*Re(y) << endl;
    cout << "Im(x)*Im(y): " << Im(x)*Im(y) << endl;
    cout << "Im(x)*Re(y): " << Im(x)*Re(y) << endl;
    cout << "Re(x)*Im(y): " << Re(x)*Im(y) << endl;
    cout << "Re(x)*Re(y)+Im(x)*Im(y): " << Re(x)*Re(y)+Im(x)*Im(y) << endl;
    cout << "Im(x)*Re(y)-Re(x)*Im(y): " << Im(x)*Re(y)-Re(x)*Im(y) << endl;
    cout << "Real part: " << (Re(x)*Re(y) + Im(x)*Im(y))/den << endl;
    cout << "Imag part: " << (Im(x)*Re(y) - Re(x)*Im(y))/den << endl;
    cout << "-------- " << endl;
*/
  return INTERVAL_COMPLEX ( (Re(x)*Re(y) + Im(x)*Im(y))/den,
													  (Im(x)*Re(y) - Re(x)*Im(y))/den );

  }
  else {
  INTERVAL den = Sqr(Re(y)) + Sqr(Im(y));
  return INTERVAL_COMPLEX ( (Re(x)*Re(y) + Im(x)*Im(y))/den,
													  (Im(x)*Re(y) - Re(x)*Im(y))/den );

  }
}
#endif // NAIVE

//-----------------------------------------------------------------
#ifdef ROKNE  // For INTERVAL_COMPLEX denominators

// Div1q is a routine used when the divisor is INTERVAL_COMPLEX 
VOID Div1q(REAL& X1, REAL& X2, REAL& Y1, REAL& Y2,
					 REAL& X1D, REAL& X2D, REAL& Y1D, REAL& Y2D){

#ifdef DEBUG
  cout << "Entered Div1q: " << endl;
  cout << "X1 " << X1 << endl;
  cout << "X2 " << X2 << endl;
  cout << "Y1 " << Y1 << endl;
  cout << "Y2 " << Y2 << endl;
#endif

  INTERVAL AU12 = Sqr(Hull(X1));
  INTERVAL AU34 = Sqr(Hull(Y1));
  INTERVAL AU56 = AU12 + AU34;

  INTERVAL AK1112 = X1/AU56;
  INTERVAL AL1112 = Y1/AU56;
  INTERVAL AU78 = Sqr(Hull(X2));  
  INTERVAL AU910 = Sqr(Hull(Y2));
  INTERVAL AU1112 = AU78+AU34;
  INTERVAL AK2122 = X2/AU1112;
  INTERVAL AL2122 = Y1/AU1112;
  AU56 = AU78 + AU910;
  INTERVAL AK3132 = X2/AU56;
  INTERVAL AL3132 = Y2/AU56;
  AU56 = AU12 + AU910;
  INTERVAL AK4142 = X1/AU56;
  INTERVAL AL4142 = Y2/AU56;
#ifdef DEBUG
  cout << AK1112 << endl;
  cout << AK2122 << endl;
  cout << AK3132 << endl;
  cout << AK4142 << endl;
  cout << AL1112 << endl;
  cout << AL2122 << endl;
  cout << AL3132 << endl;
  cout << AL4142 << endl;
#endif
  REAL A1 = Min(Min(Inf(AK1112),Inf(AK2122)),Min(Inf(AK3132),Inf(AK4142)));
  REAL A2 = Max(Max(Sup(AK1112),Sup(AK2122)),Max(Sup(AK3132),Sup(AK4142)));
  REAL B1 = Min(Min(Inf(AL1112),Inf(AL2122)),Min(Inf(AL3132),Inf(AL4142)));
  REAL B2 = Max(Max(Sup(AL1112),Sup(AL2122)),Max(Sup(AL3132),Sup(AL4142)));
#ifdef DEBUG
  cout << endl;
  cout << "A1 " << A1 << endl;
  cout << "A2 " << A2 << endl;
  cout << "B1 " << B1 << endl;
  cout << "B2 " << B2 << endl;
#endif

  X1D = A1;
  if ( (Y1 > X2) || ((X1>Y2 && Y1>0)) ){
    X2D = A2;
    Y1D = -B2;
    Y2D = -B1;
    return;
	}
  else if ( (Y2 >= X2 && X2 >= Y1 && Y1 >= X1) || (X2>Y2 && Y1>=X1)) {
		X2D = Max(A2, (1.0/Y1)/2);
		Y1D = -B2;
		Y2D = -B1;
    return;
	}
  REAL P1 = 1.0 / X1;
  REAL Q1 = P1 / 2.0;
  if (Y1 <= 0) {
	}
	else {
    X2D = A2;
    Y1D = Min(-B2,-Q1);
    Y2D = -B1;
    return;
	}
  X2D = P1;
  if (Abs(Y1)>=X1){
    Y2D = Q1;
	}
  else Y2D = -B1;
  if (Abs(Y2)>=X1) Y1D = -Q1;
  else Y1D = -B2;
  return;
}

// One-argument divisor with INTERVAL_COMPLEX denominator 
INTERVAL_COMPLEX & INTERVAL_COMPLEX::operator /= (INTERVAL_COMPLEX & x){
// /* fixme */
}

// Two-argument division with INTERVAL_COMPLEX denominator
INTERVAL_COMPLEX operator/ (REAL x, const INTERVAL_COMPLEX& y){
  REAL LP, NP, MP, PP;
  REAL L = Inf(Re(y));
  REAL M = Inf(Im(y));
  REAL N = Sup(Re(y));
  REAL P = Sup(Im(y));

  if ( (M > 0 && L >= 0) || (P > 0 && L > 0) ) 
		Div1q(L,N,M,P,LP,NP,MP,PP);
  else if ( (M > 0 && L < 0.0) || (M >= 0 && N < 0) ) {
		Div1q(M,P,-N,-L,PP,MP,LP,NP);
		MP = -MP;
    PP = -PP;
	}
  else if ( (N<0 && M < 0) || (N<=0 && P < 0) ) {
		Div1q(-N,-L,-P,-M,NP,LP,PP,MP);
		LP = -LP;
		MP = -MP;
		NP = -NP;
		PP = -PP;
	}
  else {
    Div1q(-P,-M,L,N,MP,PP,NP,LP);
		NP = -NP;
		LP = -LP;
	}
  // else ERROR: Division by zero.
  INTERVAL DReal(LP,NP), DImag(MP, PP);
  INTERVAL_COMPLEX d(DReal, DImag);
  INTERVAL_COMPLEX Result = x * d;
  return (Result);
}
INTERVAL_COMPLEX operator/ (INTERVAL x, const INTERVAL_COMPLEX& y){
  REAL LP, NP, MP, PP;
  REAL L = Inf(Re(y));
  REAL M = Inf(Im(y));
  REAL N = Sup(Re(y));
  REAL P = Sup(Im(y));

  if ( (M > 0 && L >= 0) || (P > 0 && L > 0) ) 
		Div1q(L,N,M,P,LP,NP,MP,PP);
  else if ( (M > 0 && L < 0.0) || (M >= 0 && N < 0) ) {
		Div1q(M,P,-N,-L,PP,MP,LP,NP);
		MP = -MP;
    PP = -PP;
	}
  else if ( (N<0 && M < 0) || (N<=0 && P < 0) ) {
		Div1q(-N,-L,-P,-M,NP,LP,PP,MP);
		LP = -LP;
		MP = -MP;
		NP = -NP;
		PP = -PP;
	}
  else {
    Div1q(-P,-M,L,N,MP,PP,NP,LP);
		NP = -NP;
		LP = -LP;
	}
  // else ERROR: Division by zero.
  INTERVAL DReal(LP,NP), DImag(MP, PP);
  INTERVAL_COMPLEX d(DReal, DImag);
  INTERVAL_COMPLEX Result = x * d;
  return (Result);
}
INTERVAL_COMPLEX operator/ (COMPLEX& x, const INTERVAL_COMPLEX& y){
  REAL LP, NP, MP, PP;
  REAL L = Inf(Re(y));
  REAL M = Inf(Im(y));
  REAL N = Sup(Re(y));
  REAL P = Sup(Im(y));

  if ( (M > 0 && L >= 0) || (P > 0 && L > 0) ) 
		Div1q(L,N,M,P,LP,NP,MP,PP);
  else if ( (M > 0 && L < 0.0) || (M >= 0 && N < 0) ) {
		Div1q(M,P,-N,-L,PP,MP,LP,NP);
		MP = -MP;
    PP = -PP;
	}
  else if ( (N<0 && M < 0) || (N<=0 && P < 0) ) {
		Div1q(-N,-L,-P,-M,NP,LP,PP,MP);
		LP = -LP;
		MP = -MP;
		NP = -NP;
		PP = -PP;
	}
  else {
    Div1q(-P,-M,L,N,MP,PP,NP,LP);
		NP = -NP;
		LP = -LP;
	}
  // else ERROR: Division by zero.
  INTERVAL DReal(LP,NP), DImag(MP, PP);
  INTERVAL_COMPLEX d(DReal, DImag);
  INTERVAL_COMPLEX Result = x * d;
  return (Result);
}
INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, 
														const INTERVAL_COMPLEX& y)
{
  REAL LP, NP, MP, PP;
  REAL L = Inf(Re(y));
  REAL M = Inf(Im(y));
  REAL N = Sup(Re(y));
  REAL P = Sup(Im(y));

  if ( (M > 0 && L >= 0) || (P > 0 && L > 0) ) 
		Div1q(L,N,M,P,LP,NP,MP,PP);
  else if ( (M > 0 && L < 0.0) || (M >= 0 && N < 0) ) {
		Div1q(M,P,-N,-L,PP,MP,LP,NP);
		MP = -MP;
    PP = -PP;
	}
  else if ( (N<0 && M < 0) || (N<=0 && P < 0) ) {
		Div1q(-N,-L,-P,-M,NP,LP,PP,MP);
		LP = -LP;
		MP = -MP;
		NP = -NP;
		PP = -PP;
	}
  else {
    Div1q(-P,-M,L,N,MP,PP,NP,LP);
		NP = -NP;
		LP = -LP;
	}
  // else ERROR: Division by zero.
  INTERVAL DReal(LP,NP), DImag(MP, PP);
  INTERVAL_COMPLEX d(DReal, DImag);
  INTERVAL_COMPLEX Result = x * d;
  return (Result);
}
#endif
//--------------------------------------------------------------

/*

// This code is modelled after PROFIL, but it isn't working as well as
// the naive approach!!!

INTERVAL_COMPLEX operator/ (const INTERVAL_COMPLEX& x, 
														const INTERVAL_COMPLEX& y)
{
//  a+bi   (ac + bd) + (bc - ad)i
//  ---- = ----------------------
//  c+di        c*c + d*d

  // I should check this!!!!

  INTERVAL r,s,den;

  if (0 <= Re(y)  && !(0 <= Im(y))) {
    r = Re(y) /Im(y);
    den = Im(y) + r * Re(y);
    return INTERVAL_COMPLEX( (Re(x)*r + Im(x))/den, (Im(x)*r - Re(x))/den ); 
	}
  else if (0 <= Im(y) && !(0 <= Re(y))){
    r = Im(y)/Re(y);
    den = Re(y) + r * Im(y);
    return INTERVAL_COMPLEX( (Re(x) + r*Im(x))/den, (Im(x) - r*Re(x))/den );
	}
  else {
    r = Re(y) /Im(y);
		s = Im(y) /Re(y);
    if ( Sup(Abs(r)) < Sup(Abs(s)) ){
      den = Im(y) + r * Re(y);
      return INTERVAL_COMPLEX( (Re(x)*r + Im(x))/den, (Im(x)*r - Re(x))/den ); 
		}
		else {
      den = Re(y) + r * Im(y);
      return INTERVAL_COMPLEX( (Re(x) + r*Im(x))/den, (Im(x) - r*Re(x))/den );
		}
	}

}
*/

INTERVAL Abs (const INTERVAL_COMPLEX& x)
{
/*
  1 2 3
    |
  4-5-6
    |
  7-8-9
*/
  INTERVAL Result;

  // Boxes 1, 2 and 3
  if (Inf(Im(x))>=0.0){
    // Box 1
    if (Sup(Re(x))<=0.0){
			Result = Hull( Sqr(Inf(Re(x)))+Sqr(Sup(Im(x))),
               Sqr(Sup(Re(x)))+Sqr(Inf(Im(x))));
    }
    else if (Inf(Re(x))>=0.0){ 
      // Box 3
			Result = Hull( Sqr(Sup(Re(x)))+Sqr(Sup(Im(x))),
			          Sqr(Inf(Re(x)))+Sqr(Inf(Im(x))));
		}
    else {
			// Box 2
			Result = Hull( Sqr(Max(Sup(Re(x)),Abs(Inf(Re(x)))))+Sqr(Sup(Im(x))),
                 0.0+Sqr(Inf(Im(x))));
 		}
	}
  // Box 7,8 or 9
  else if (Sup(Im(x))<=0.0){  
    if (Sup(Re(x))<=0.0){  
      // Box 7
      Result = Hull( Sqr(Inf(Re(x)))+Sqr(Inf(Im(x))),
                 Sqr(Sup(Re(x)))+Sqr(Sup(Im(x))));
		}
		else if (Inf(Re(x))>= 0.0){  
      // Box 9
			Result = Hull( Sqr(Sup(Re(x)))+Sqr(Inf(Im(x))),
			          Sqr(Inf(Re(x)))+Sqr(Sup(Im(x))));
		}
		else {  
      // Box 8
			Result = Hull( Sqr(Sup(Im(x))),
			           Sqr(Max(Abs(Inf(Re(x))),Sup(Re(x))))+Sqr(Inf(Im(x))));
		}
	}
  // Box 4,5 or 6
  else {  
    if (Sup(Re(x))<=0.0) { 
      // Box 4
			Result = Hull( Sqr(Max(Abs(Inf(Im(x))),Sup(Im(x))))+Sqr(Inf(Re(x))),
			          Sqr(Sup(Re(x)))+0.0);
		}
		else if (Inf(Re(x))>=0.0) {  
      // Box 6
			Result = Hull( Sqr(Max(Abs(Inf(Im(x))),Sup(Im(x))))+Sqr(Sup(Re(x))),
                 Sqr(Inf(Re(x)))+0.0);
		}
		else { 
      // Box 5
      Result = Hull( Sqr(Max(Abs(Inf(Re(x))),Sup(Re(x))))
                +Sqr(Max(Abs(Inf(Im(x))),Sup(Im(x)))),
              0.0);
		}
	}

  return Sqrt(Result);
}

ostream & operator << (ostream & os, const INTERVAL_COMPLEX& x)
{
  return os << '(' << Re (x) << ',' << Im (x) << ')';
}

// The following code is only used to force Constants.C to be
// always included in the executable code.

extern VOID RegisterConstants ();

//VOID RegisterComplex () { RegisterConstants (); }
