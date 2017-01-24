// bivf.cc
// Studentized chi-square distributions


#include "MiscFunctions.h"    // For the Random() and Rand01() functions.
#include "Functions.h"      // For Cos, Sqr, Constant::Pi
#include "UnconstrainedOpt.h"
#include "Interval.h"         // for Intersection
#include "Constants.h"        // for Machine::PosInfinity
#include "TestMatrices.h"     // For "matrix" stuff.
#include "Utilities.h"
#include "IntegerVector.h"


// Prototypes
//INTERVAL Gamma(INTERVAL x);
//INTERVAL IncGammaCDF(REAL Alpha, INTERVAL x);
//INTERVAL IncGammaCDF(INTERVAL Alpha, INTERVAL x);
//INTERVAL LogGamma(INTERVAL x);

// Functions definitions

INTERVAL LogGamma(INTERVAL x){
  // Calculate the Log(Gamma(x))
  // This function has been checked for several dozen values and appears OK.

  INTERVAL Result;
  INT k;

  if (x == Hull(0.5))
    Result = Log(Hull(1.772453850905516));
  else if (Sup(x) < 10){
    static INTERVAL_VECTOR ck(26);
    ck(1) = Hull(1.0);
    ck(2) = Hull(0.57721566490153,0.57721566490154);
    ck(3) = Hull(-.65587807152026,-.65587807152025);
    ck(4) = Hull(-.042002635034096,-.042002635034095);
    ck(5) = Hull(.16653861138229,.16653861138230);
    ck(6) = Hull(-.0421977345555444,-.0421977345555442);
    ck(7) = Hull(-.0096219715278771,-.0096219715278769);
    ck(8) = Hull(.0072189432466629,.0072189432466631);
    ck(9) = Hull(-.0011651675918592,-.0011651675918590);
    ck(10) = Hull(-.000215241674115,-.0002152416741149);
    ck(11) = Hull(.0001280502823881,.0001280502823883);
    ck(12) = Hull(-.201348547808e-4,-.201348547806e-4);
    ck(13) = Hull(-.12504934822e-5,-.12504934820e-5);
    ck(14) = Hull(.11330272319e-5,.11330272321e-5);
    ck(15) = Hull(-.2056338418e-6,-.2056338417e-6);
    ck(16) = Hull(.61160949e-8,.61160951e-8);
    ck(17) = Hull(.50020074e-8,.50020076e-8);
    ck(18) = Hull(-.11812747e-8,-.11812745e-8);
    ck(19) = Hull(.1043426e-9,.1043428e-9);
    ck(20) = Hull(.77822e-11,.77824e-11);
    ck(21) = Hull(-.36969e-11,-.36967e-11);
    ck(22) = Hull(.5099e-12,.5101e-12);
    ck(23) = Hull(-.207e-13,-.205e-13);
    ck(24) = Hull(-.55e-14,-.53e-14);
    ck(25) = Hull(.13e-14,.15e-14);
    ck(26) = Hull(.09e-15,.20e-15);

    // First we use Gamma(x) = (x-1) * Gamma(x-1) repeatedly until 1 < x <= 2
    // Note that are summing logs, since we return Log(Gamma(x))
    INTERVAL Sum = 0.0;
    while (Sup(x)>2.0) {
      x -= 1.0;
      Sum += Log(x);
    }    

    // Then use Abromowitz/Stegun, 6.1.34, p. 256, i.e.
    // 1/Gamma(x) = sum_1^infinity c_k * x^k
    INTERVAL TaylorSum = 0.0;
    for (k=1; k<26; k++)
      TaylorSum += ck(k) * Power(x,k) ;

    // After summing the first 25 terms of the series, we calculate the
    // truncation error, then add & subtract it from the Taylor series sum 
    INTERVAL r = ( ck(26)*Power(x,26) ) / ( ck(25)*Power(x,25) ) ;
    INTERVAL ErrorTerm = ck(26)*Power(x,26) / (1-r);
    TaylorSum += Hull(-ErrorTerm,ErrorTerm);
    Result = Sum + Log(1/TaylorSum) ;    
  }

  else{ // Sup(x) >= 10
    static INTERVAL_VECTOR Bern(30);
    Bern(1) = Hull(0.16666666666666,0.16666666666667);
    Bern(2) = Hull(-.033333333333334,-.033333333333333);
    Bern(3) = Hull(.02380952380952380,.02380952380952381);
    Bern(4) = Hull(-.033333333333334,-.033333333333333);
    Bern(5) = Hull(.075757575757575,.075757575757576);
    Bern(6) = Hull(-.25311355311354,-.25311355311355);
    Bern(7) = Hull(1.166666666666665,1.166666666666666);
    Bern(8) = Hull(-7.092156862745096,-7.092156862745097);
    Bern(9) = Hull(54.97117794486214,54.97117794486215);
    Bern(10) = Hull(-529.12424242424243,-529.12424242424242);
    Bern(11) = Hull(6192.123188415797,6192.123188415798);
    Bern(12) = Hull(-8650.253112553118,-8650.253112553117);
    Bern(13) = Hull(1425517.1666666667,1425517.1666666668);
    Bern(14) = Hull(-27298231.167816095,-27298231.167816094);
    Bern(15) = Hull(601580873.9006424,601580873.9006425);
    Bern(16) = Hull(-15116315767.092158,-15116315767.092157);
    Bern(17) = Hull(429614643061.16669,429614643061.16670);
    Bern(18) = Hull(-13711655205088.335,-13711655205088.334);
    Bern(19) = Hull(488332318973593.19,488332318973593.20);
    Bern(20) = Hull(-19296579341940068.,-19296579341940067.);
    Bern(21) = Hull(841693047573682690.,841693047573682691.);
    Bern(22) = Hull(-40338071854059454e3,-40338071854059453e3);
    Bern(23) = Hull(21150748638081993e5,21150748638081994e5);
    Bern(24) = Hull(-12086626522296526e7,-12086626522296525e7);
    Bern(25) = Hull(75008667460769642e8,75008667460769643e8);
    Bern(26) = Hull(-50387781014810682e10,-50387781014810681e10);
    Bern(27) = Hull(36528776484818122e12,36528776484818123e12);
    Bern(28) = Hull(-28498769302450882e14,-28498769302450881e14);
    Bern(29) = Hull(23865437499683627e16,23865437499683628e16);
    Bern(30) = Hull(-21399949257225338e17,-21399949257225337e17);

    // For larger values of x, use Abromowitz/Stegun, 6.1.40, p. 257, i.e.
    // LogGamma(x) = (x-.5)ln(x) - x + .5ln(2pi) 
    //             + sum_1^infinity B_2m / [ 2m*(2m-1)x^(2m-1) ]

    // First the non-series part
    INTERVAL PartA  = (x-0.5)*Log(x) - x + 0.5 * Log(Constant::TwoPi);

    // Then the series part
    INTERVAL Term = 0.0;
    INTERVAL Sum = 0.0;
    for (k=1; k<= 30; k++){
      Term = Bern(k) / (2*k*(2*k-1)*Power(x,(2*k-1)) );
      if (k < 30){
        Sum += Term;
      }
      else {   // For k=30, Term is the remainder, so add & subtract it
        Sum += Hull(-Term, Term);   // See the note in 6.1.42
      }
    }
    Result = PartA + Sum;
  }  // end of else Sup(x) >= 10

  return (Result);
}

INTERVAL Gamma(INTERVAL x){
  // Calculate Gamma(x)
  // The function has been moderately-well checked.

  return (Exp(LogGamma(x)));
}

VOID Derivs(INTERVAL x, INT k, INTERVAL a, INTERVAL_VECTOR& U, 
      INTERVAL_VECTOR& V, INTERVAL_VECTOR& Coeff){
  // In the Taylor expansion, evaluate the k^th derivative terms of exp(-x) 
  // and x^(a-1).  See Wang, p. 15

  // Note that the first vector element is the "zero" derivative
  // so care must be taken with indices.  I.e. the k^th derivative is the
  // (k+1)^th vector element

  // k^th "derivative" term of exp(-x), i.e. (-1)^k exp(-x) / k! 
  U(k+1) = U(k)/(-k);

  // k^th derivative of x^(a-1)
  V(k+1) = (a-k)*V(k) / (k * x);

  INTERVAL Coef = 0.0;
  for (INT j=1; j<=k+1; j++)
    Coef += U(j)*V(k-j+2);
  
  Coeff(k+1) = Coef;
  return;
}

INTERVAL TaylorSeries(INTERVAL X, INTERVAL C, INTERVAL Alpha, INTERVAL H){
  // Employ Corliss & Rall's algorithm for Taylor Series of Inc. Gamma fn.
  // See Wang's thesis, page 22.
  INTERVAL_VECTOR Uc(1000), Vc(1000), Ux(1000), Vx(1000);
  INTERVAL_VECTOR Fc(1000), Fx(1000);
  INTERVAL_VECTOR Hpart(1000);
  INTERVAL Sum, Newj, Oldj, jIntersect;

  // First initialize the "zero" derivative of Exp(-C), C = midpoint of C
  Uc(1) = 1/Exp(C);
  // and the "zero" derivative of C^(Alpha-1)
  Vc(1) = Power(C, Alpha-1);

// I don't think this "if" is ever satisfied...
  // Now initialize the "zero" derivatives of Exp(-X) and X^(Alpha-1)
  // If the lower endpoint of X is 0, we move the lower endpoint to be the
  // same as the upper endpoint and thereby avoid division by zero in 
  // calculating the next Vx().
  if (Inf(X) == 0){
    X = Hull(Sup(X));
    Ux(1) = Hull(1/Exp(X), 1.0);
    Vx(1) = Hull(0.0, Power(X, Alpha-1));
  }
  else {
    Ux(1) = 1/Exp(X);
    Vx(1) = Power(X,Alpha-1);
  }

  // Now calculate the "zero" derivative coefficient in the taylor series
  Fc(1) = Uc(1) * Vc(1);
  Fx(1) = Ux(1) * Vx(1);

  INTERVAL hpow = H;
  Hpart(1) = hpow/1.0;      // Hpart(i) = H^i/i!
  hpow *= H;
  Hpart(2) = hpow/2.0;

  // Determine J_n for n=1 (and so i=0)
  INT n=1;  // Not really needed, just keeping things straight
  // First we need to get the next Fx().  While we're at it, update Fc().
  Derivs(C,1,Alpha,Uc,Vc,Fc); 
  Derivs(X,1,Alpha,Ux,Vx,Fx);
  Sum = Fc(1) * Hpart(1);  // For n=1, i=0. 
  INTERVAL R1 = Fx(n+1) * Hpart(n+1);      // Fx(2) is first derivative
  Oldj = 2.0 * Sum + (R1-R1);

  // Compute successive J_n's until the interval is narrow enough
  // or has ceased decreasing
  BOOL Converge = 0;

  while (!Converge){
    n++;

    // Now figure the next j.  First update the Fc() and Fx() coefficients
    Derivs(C,n,Alpha,Uc,Vc,Fc);
    Derivs(X,n,Alpha,Ux,Vx,Fx);
    hpow *= H;
    Hpart(n+1) = hpow/(float)(n+1);

    // Compute the remainder term
    if ( n/2.0 == (int)(n/2.0) ){  // even n
      // Note we use Fx(n+1) to get the n^th derivative coefficient!
      Newj = 2.0 * (Sum + Fx(n+1) * Hpart(n+1) );
    }
    else{  // odd n
      // Note we use Fc(n) to get the (n+1)-1 ^th coeffecient!
      Sum += Fc(n) * Hpart(n);  // Update Sum for odd values of n
      R1 = Fx(n+1) * Hpart(n+1);
      Newj = 2.0 * Sum + (R1-R1);
    }
    if (Newj == Oldj)  Converge = 1;

    // Jint = Newj (intersect) Oldj;
    Intersection(jIntersect,Newj,Oldj);

    Oldj = jIntersect;

  }

  // Now divide by Gamma(Alpha) and return
  return (jIntersect/Gamma(Alpha));
}

INTERVAL IncGammaCDF(REAL Alpha, INTERVAL x){
// Incomplete Gamma function
  // Note, although the declaration says Alpha is a REAL, it is really
  // required to be a machine floating number.  This function is
  // called to evaluate IncGammaCDF when Alpha is one of the endpoints
  // of an interval.
  // Although REAL Alpha will always be promoted to INTERVAL Alpha,
  // a separate variable is declared for emphasis.  Since Alpha is a machine
  // number, Hull(Alpha) will always be degenerate width.

  INTERVAL Result;

  if (Sup(x) <= 0)
    Result = Hull(0.0);

  else if (Alpha == (int)Alpha){  // Case 1.  Alpha is an integer ----------

    INT i;
    INTERVAL Term = 1.0;  // 0^th Term of summation
    INTERVAL Comp = 1.0; // Complement of result

    for (i=1; i<Alpha; i++){
      Term *= x/i;
      Comp += Term;
    }
    Comp /= Exp(x);
    Result = 1 - Comp;
  }

  
  else if (Inf(x) > 1.5){ // Case 2.  Alpha not integer, x > 1.5 ------------

    INT n = (int)Alpha;  // n = integral part of alpha
    INTERVAL b = Alpha - n;

    // Priming information
    INTERVAL p_kminus2 = 1.0, p_kminus1 = 0.0;
    INTERVAL q_kminus2 = 0.0, q_kminus1 = 1.0;
    INTERVAL c_kminus1 = 0;
    INTERVAL a_k, b_k, c_k, p_k, q_k;

    // The first convergent
    INT k = 1;  // Some c.f's start at k=0, but k=1 is more convenient here
    a_k = 1.0;
    b_k = 1.0;
    p_k = b_k * p_kminus1 + a_k * p_kminus2;
    q_k = b_k * q_kminus1 + a_k * q_kminus2;
    c_k = p_k / q_k;

    BOOL Converge = 0;  
    while (!Converge) {
      k++;

      if ( k/2.0 == (int)(k/2) ) {  // Even k

        a_k = (int)(k/2) - b;  // Note that we're using integer division
        b_k = x;
      }
      else {    // Odd k
        a_k = (int)(k/2);
        b_k = 1;
      }

      p_kminus2 = p_kminus1;    
      p_kminus1 = p_k;
      q_kminus2 = q_kminus1;
      q_kminus1 = q_k;
      c_kminus1 = c_k;
      p_k = b_k * p_kminus1 + a_k * p_kminus2;
      q_k = b_k * q_kminus1 + a_k * q_kminus2;
      // Rescale the problem from time to time.  Prevents overflow.
      if (Sup(p_k) > 1.0e120){
        p_k /= 1.0e100; p_kminus1 /= 1.0e100;  p_kminus2 /= 1.0e100;
        q_k /= 1.0e100; q_kminus1 /= 1.0e100;  q_kminus2 /= 1.0e100;
      }
      c_k = p_k / q_k;
    
      Converge = (c_k == c_kminus1);
    }
    INTERVAL ContFrac = c_k;

    // Now for the rational fractional part.  This is calculated if n>0.
    INTERVAL RatFrac = 0.0;
    if (n > 0) {
      INT i = 0;
      INTERVAL Term = 1.0/b;
      RatFrac += Term; 
      for (i = 1; i<n; i++){
        Term *= x/(b+i);
        RatFrac += Term;
      }
    }

    // Now add the two pieces together and multiply by the common constant
    INTERVAL Complement = (Power(x,b-1)*ContFrac + Power(x,b)*RatFrac);
    Complement /= Exp( x + LogGamma(b) );
    Result = 1 - Complement;
  }

  else {  // Case 3.  Inf(x) <= 1.5 ---------------------------------
    // Split the integgral into two pieces.  Use an alternating series for the
    // left integral and Taylor Series for the right integral.
    // This split is necessary to avoid the singularity in the derivatives
    // of t^(a-1)

    // Determine midpoint of 0 to X
    INTERVAL MidPt = x/2.0;
    INTERVAL Md = (MidPt + x) / 2.0;
    // Determine interval hull
    INTERVAL y = Md - MidPt;
    INTERVAL z = x - Md;
    INTERVAL h = Hull( Hull(Inf(y),Sup(z)), Hull(Sup(y),Sup(z)) );
    INTERVAL x0 = x/2.0; 
    
    // First we compute 1/Gamma(a) * int_0^x0 e^-t t^(a-1)
    // This is where the call to "modsing" is made
    
    // EvenSum is the sum of the series through an even term
    INT i = 0;
    INTERVAL Term = 1;   // Term is (-x0)^i / (i!)
    INTERVAL EvenSum = Term/(Alpha+i);       
    
    // OddSum is the sum of the series through an odd term 
    i = 1;
    Term *= -x0/i;       
    INTERVAL OddSum = EvenSum + Term/(Alpha+i);  // First term i=1
    REAL PrevWidth = 100.0;
    REAL Width = 1.0;
    INT MaxIter = 1000;
    
    // Since this is an alternating series, we calculate two terms at
    // a time and use these as the upper and lower bounds on the sum
    // until the width of the interval doesn't change.
    while ((Width < PrevWidth) && i < MaxIter){             // F
      i++;  // After this line, i is even, so calculate an even term.
      Term *= -x0/i;
      EvenSum = OddSum + Term/(Alpha + i);
      
      i++;  // After this line, i is odd
      Term *= -x0/i;
      OddSum = EvenSum + Term/(Alpha + i);
      
      PrevWidth = Width;
      Width = Sup(EvenSum) - Inf(OddSum);
    } 
    INTERVAL AltSeries = Hull(Inf(OddSum), Sup(EvenSum));
    // Now premultiply by:  x0^a / Gamma(a)
    AltSeries *= Exp( Alpha*Log(x0) - LogGamma(Alpha) );
    
    // OK, now we compute 1/Gamma(a) * int_x^x0 e^-t t^(a-1) by using
    // a Taylor series.   Left term of RHS of (5.12)
    // First thing we do is pull a trick and change the lower endpoint
    // of the interval 
    x = Hull(Inf(MidPt), Sup(x));
    INTERVAL Taylor = TaylorSeries(x, Md, Alpha, h);

    // Lastly, add the two pieces together
    Result = AltSeries + Taylor;
      
  }  // End of if alpha < 1.0

  // Before returning the value of the incomplete Gamma function
  // we make sure it is contained in [0,1].
  if (Inf(Result) < 0) Result = Hull(0, Sup(Result));
  if (Sup(Result) > 1) Result = Hull(Inf(Result), 1);
  return (Result);
}

INTERVAL IncGammaCDF(INTERVAL Alpha, INTERVAL X){
  // Incomplete Gamma function
  // Since this is a CDF, it is monotone increasing, and it is sufficient
  // to evaluate the function at each endpoint and take the hull of the
  // result.  If Alpha is degenerate, then we only need to call the 
  // scalar function once.

  INTERVAL Result = 0.0;
  if (Inf(Alpha) == Sup(Alpha))
    Result = IncGammaCDF(Inf(Alpha), X);
  else
    Result = Hull( IncGammaCDF(Inf(Alpha),X),IncGammaCDF(Sup(Alpha),X) );
  return (Result);
}

/*
INTERVAL f(REAL x, REAL df){
  INTERVAL X = Hull(x);
  INTERVAL R;
  R = IncGammaCDF(df, X);
  return (R);
}
INTERVAL fp(INTERVAL X, REAL df){
  INTERVAL R;
   INTERVAL Df = Hull(df);
  R = Exp( -X + (Df-1)*Log(X)  - LogGamma(Df) );
  return (R);
}
*/ 
/*
INTERVAL Zero(INTERVAL X, REAL df, REAL pval){
  // Needs f(REAL x) & fp(INTERVAL X)
  
  cout << "Input pval: " ;
  cin >> pval;

  BOOL Done = 0;
  INTERVAL Xn = X, Xnp1, Xtemp;
  REAL xm;

  if (0 <= fp(Xn, df)) Done = 1;
  if (Sup( (f(Inf(Xn),df)-pval)*(f(Sup(Xn),df)-pval) )>= 0) Done = 1;

  while (!Done){
    
    // Now compute the next iterate
    xm = (Inf(Xn)+Sup(Xn)) / 2.0;
    Xtemp = xm - (f(xm, df)-pval) / fp(Xn, df);

    Intersection(Xnp1, Xn, Xtemp);

    if (Xn == Xnp1) Done = 1;

    // Prepare for next iteration                  
    Xn = Xnp1;
    if (0 <= fp(Xn, df)) Done = 1;
    if (Sup( (f(Inf(Xn),df)-pval)*(f(Sup(Xn),df)-pval) )>= 0) Done = 1;
  }

  return (Xn);
}
*/

INTERVAL g(INTERVAL X, INT m, INT n, INTERVAL d1, INTERVAL d2, 
            INTERVAL Rho, INT j){
  // This is the integrand of interest
  INTERVAL Result = 0.0;
  INTERVAL HalfN = Hull(n/2.0);
  INTERVAL Rho2 = Sqr(Rho);

  if (d1 == d2){
    Result = Exp(-X/2.0)*Power(X, HalfN-1) * 
               Sqr(IncGammaCDF(m/2.0+j, d1*m*X/(2.0*n*(1-Rho2))));
  }
  else{
    Result = Exp(-X/2.0)*Power(X, HalfN-1) * 
      IncGammaCDF(m/2.0+j, d1*m*X/(2.0*n*(1-Rho2))) *
      IncGammaCDF(m/2.0+j, d2*m*X/(2.0*n*(1-Rho2)));
  }
  Result /= Power(2.0, n/2.0)*Gamma(HalfN);
  return (Result);  
}

INTERVAL gpp(INTERVAL X, INT m, INT n, INTERVAL d1, INTERVAL d2, INTERVAL Rho, INT j){
  // Second derivative of g
  INTERVAL Result = 0.0;
  INTERVAL HalfN = Hull(n/2.0);
  INTERVAL C = m/2.0 +j;
  INTERVAL Rho2 = Sqr(Rho);
  INTERVAL K1 = d1*m / ( 2*n*(1-Rho2) ), K2 = d2*m / ( 2*n*(1-Rho2) );
  INTERVAL G1 = IncGammaCDF(C,K1*X);
  INTERVAL G2 = IncGammaCDF(C,K2*X);
  INTERVAL GC = Gamma(C);

  Result = G1 * G2 * Exp(-X/2.0) * Power(X,HalfN-3) * 
             (Sqr(X)/4.0 -(HalfN-1)*X +(HalfN-1)*(HalfN-2));
  Result += G1 * Exp(-(0.5+K2)*X) * Power(X,HalfN+C-3) *
             (-X -K2*X +n+C-3) * Power(K2,C) / GC;
  Result += G2 * Exp(-(0.5+K1)*X) * Power(X,HalfN+C-3) *
             (-X -K1*X +n+C-3) * Power(K1,C) / GC;
  Result += 2.0 * Power(K1*K2,C) * Power(X,HalfN+C+C-3) * 
             Exp(-(0.5+K1+K2)*X) / Sqr(GC);
 
  Result /= Power(2.0, n/2.0)*Gamma(HalfN);  
  return (Result);  
}

INTERVAL NewtonCotes(INTERVAL X, INT m, INT n, 
             INTERVAL d1, INTERVAL d2, INTERVAL Rho, INT j, INT count){
  // Calculate the Newton-Cotes approximation to int_a^b g dx
  // and include the error term
  INTERVAL Result = 0;
  REAL a = Inf(X), b = Sup(X);
  INTERVAL H = (b-a)/count;
  INTERVAL Point_i = Hull(a);
  
  INT i = 0;
  Result = g(Point_i, m, n, d1, d2, Rho, j) / 2.0;

  for(i=1; i<count; i++){
    Point_i = a + H*i;
    Result += g(Point_i, m, n, d1, d2, Rho, j);    
//    cout << "g" << i << " " << g(Point_i, m, n, d1, d2, Rho, j) << endl;
//    cout << "Result: " << Result << endl;
  }

  i = count;
  Point_i = Hull(b);
  Result += g(Point_i, m, n, d1, d2, Rho, j) / 2.0;

  Result *= H;

  Result += -Power(b-a,3)*gpp(X, m, n, d1, d2, Rho, j)/(12*Sqr(count));

  return (Result);
}

INTERVAL B(INT m, INT n, INTERVAL d1, INTERVAL d2, INT j, INTERVAL Rho,
           REAL Eps1, REAL Eps2, INT nccount, REAL Width){
  INTERVAL Result, X_k;
  INTERVAL NC, LeftProb, RightProb, MidProb;
  INT k;

  // First figure the enclosure for the lower tail probability
  /* FIXME */  
  // Bound can be improved
  LeftProb = Hull(0.0, Eps1);

  // Then the enclosure for the upper tail probability
  RightProb = Hull(0.0, 1-IncGammaCDF(n/2.0, Eps2/2.0));


  // then the enclosure for the middle area.
  MidProb = Hull(0.0);
  REAL L = Eps1;
  REAL R = Eps1 + Width;
  if (R >= Eps2) {
    R = Eps2;
    X_k = Hull(L,R);
//    cout << "X_1" << "  " << X_k << endl;
    NC = NewtonCotes(X_k, m, n, d1, d2, Rho, j, nccount);    
//    cout << "NC: " << NC << endl;
    MidProb += NC; 
  }
  else {
    k = 0;
    while (R < Eps2) {
      k++;
      X_k = Hull(L,R);
//      cout << "X_" << k << "  " << X_k << endl;
      NC = NewtonCotes(X_k, m, n, d1, d2, Rho, j, nccount);    
//      cout << "NC: " << NC << endl;
      MidProb += NC; 
      L = R;             // Now move to the next interval
      R += Width;
    }
    R = Eps2;  // Move R back a bit to catch the last part
    k++;
    X_k = Hull(L,R);
//    cout << "X_" << k << "  " << X_k << endl;
    NC = NewtonCotes(X_k, m, n, d1, d2, Rho, j, nccount);    
//    cout << "NC: " << NC << endl;
    MidProb += NC; 
  }

  Result = LeftProb + MidProb + RightProb;
  cout << "LeftProb: " << LeftProb << endl;
  cout << "MidProb: " << MidProb << endl;
  cout << "RightProb: " << RightProb << endl;
  cout << "B_" << j << " : " << Result << endl;
  return Result;
}

INTERVAL BivF(INT m, INT n, INTERVAL d1, INTERVAL d2, INTERVAL Rho,
              REAL Eps1, REAL Eps2, INT nccount, REAL Width){
  // Calculate the distribution of the bivariate f variate
  // See Krishnaiah, Handbook of Statistics
  INTERVAL B_j, ProbSum = 0.0;
  INTERVAL ErrBoundSum = 0.0, ErrBound = 0.0, PrevErrBound = 0.0, Term;
  INTERVAL Rho2 = Sqr(Rho);
  INTERVAL HalfM = Hull(m/2.0);
  INTERVAL GammaHalfM = Gamma(HalfM);
  INTERVAL ErrBoundConst = Power(1-Rho2,HalfM);
  INTERVAL PartialProb = 0.0, PrevPP = 0.0;
  BOOL Converge = 0;
  INT j = 0;
  INTERVAL Numer = 1.0, Denom = 1.0;
  while (!Converge){
    B_j = B(m,n,d1,d2,j,Rho,Eps1,Eps2,nccount,Width);

    // I made substantial improvements in accuracy by eliminating the calls
    // to the Gamma function
    // Term = Gamma(m/2+j)*Rho2 / ( Gamma(m/2)*j! )
    Term = Numer / Denom;
    ProbSum += Term * B_j;
    ErrBoundSum += Term;
    ErrBound = 1-ErrBoundConst * ErrBoundSum;

    if (PrevErrBound == ErrBound) 
      Converge = 1;
    else if ( (j>0) && 
      (Abs(Inf(PrevPP)-Inf(PartialProb))<1e-10) && 
      (Abs(Sup(PrevPP)-Sup(PartialProb))<1e-10) )
      Converge = 1;
    else{
      // Now prepare for next iteration of while loop
      j++;
      // Rescale from time to time to prevent overflow
      if (Sup(Numer) > 1.0e300 || Sup(Denom) > 1.0e300) {
        Numer /= 1.0e150;
        Denom /= 1.0e150;
      }
      Numer *= (m + 2*(j-1)) * Rho2;
      Denom *= j * 2;
      PrevPP = PartialProb;
      PrevErrBound = ErrBound;
      PartialProb = ProbSum * Power(1-Rho2,HalfM);
      PartialProb = Hull(Inf(PartialProb), Sup(PartialProb)+ErrBound);
      cout << endl;
    }
  }

  // First adjust for the constant multiplier to get the probability
  INTERVAL Prob = ProbSum * Power(1-Rho2,HalfM);
  // Then enlarge the right endpoint by the error bound
  Prob = Hull(Inf(Prob), Sup(Prob)+ErrBound);

  cout << endl << "FinalProb: " << Prob << endl;
  return (Prob);
}

INTERVAL InvBivF(REAL x_iminus1, REAL x_i, INTERVAL Prob, 
         INT m, INT n, INTERVAL Rho, REAL Eps1, REAL Eps2, 
         INT nccount, REAL Width){
  INT i=0;  
  BOOL Done = 0;
  REAL Epsilon = 1e-15;
  REAL L = x_iminus1;
  REAL U = x_i;
  INTERVAL FL = BivF(m,n,L,L,Rho,Eps1,Eps2,nccount,Width) - Prob;
  INTERVAL FU = BivF(m,n,U,U,Rho,Eps1,Eps2,nccount,Width) - Prob;
  INTERVAL F_i = FU;
 
  REAL Q, PrevQ, PrevL, PrevU, x, x_iplus1;
  INTERVAL Fx, F_iplus1, IX_c, IX_i, IX_iplus1, X_iplus1;

  if (Sup(FL) < 0  &&  Inf(FU) > 0) ;
    // values are satisfactory
  else {
    cout << "==========================================" << endl;
    cout << "Initial values unsatisfactory." << endl;
    cout << "FL : " << FL << endl;
    cout << "FU : " << FU << endl;
    Done = 1;
  }

  while (!Done){
    i++;

    IX_i = Hull(L,U);
    IX_c = U - (U-L)/(1-FL/FU);
    x_iplus1 = 0.5 * (Inf(IX_c) + Sup(IX_c));
    F_iplus1 = BivF(m,n,x_iplus1,x_iplus1,Rho,Eps1,Eps2,nccount,Width) - Prob;

    if (Inf(F_iplus1) > 0) { 
      // Use the left sub-interval for the next iteration, so move the
      // right endpoint
      U = x_iplus1;
      FU = F_iplus1;
      if ( Inf(F_i * F_iplus1) > 0 )  FL /= 2.0;
      IX_iplus1 = Hull(L,U);
    }
    else if (Sup(F_iplus1) < 0) {
      // Use the right sub-interval for the next iteration, so move the
      // left endpoint
      L = x_iplus1;
      FL = F_iplus1;
      if ( Inf(F_i * F_iplus1) > 0 ) FU /= 2.0;
      IX_iplus1 = Hull(L,U);
    }
    else {
      // Don't know if X_iplus1 is the left or right sub-interval, so exit
      Done = 1;
    }

    if (0 <= (FU - FL)){
      // At the next iteration, the denominator will contain 0, so exit
      Done = 1;
    }
    // Prepare for next iteration
    F_i = F_iplus1;
  }

  cout << endl;
//  cout << "Bracket secant : " << i << " " << Hull(L,U) << endl;
 
  // Now switch to bisection search
  // First narrow the interval by moving the right endpoint
  PrevQ = PrevU = -1;
  Q = L;
  while ( (Q != PrevQ) || (U != PrevU) ) {
    PrevQ = Q;
    PrevU = U;
    i++;
    x = (Q + U) / 2.0;
    Fx = BivF(m,n,x,x,Rho,Eps1,Eps2,nccount,Width) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) Q = L = x;
    else Q = x;
  }
//  cout << "Right bisection: " << i << " " << Hull(L,U) << endl;

  // Then narrow the interval by moving the left endpoint
  PrevQ = PrevL = -1;
  Q = U;
  while ( (Q != PrevQ) || (L != PrevL) ) {
    PrevQ = Q;
    PrevL = L;
    i++;
    x = (L + Q) / 2.0;
    Fx = BivF(m,n,x,x,Rho,Eps1,Eps2,nccount,Width) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) L = x;
    else Q = x;
  }
//  cout << "Left  bisection: " << i << " " << Hull(L,U) << endl;

  return (Hull(L,U));
}



INT main(){
  cout.precision(15);
  INTERVAL X,Rho, Alpha, d1, d2;
  REAL x_iminus1, x_i, pval;
  INT j,m,n,nccount,Choice;  
  REAL Eps1, Eps2, Width;
  cout << endl;
  cout << "1  Inc Gamma CDF " << endl;
//  cout << "2  Find truncation point" << endl;
//  cout << "3  Test g gpp" << endl;
  cout << "4  BivF CDF " << endl;
  cout << "5  Newton Cotes: " << endl;
  cout << "6  B: " << endl;
  cout << "7  Find Critical Point: " << endl;
  cout << "8  Find row of Critical Points: " << endl;
  cout << endl;
  cout << "Enter choice: " ;

  cin >> Choice;

  if (Choice==1){
    cout << "Input interval X: " ;
    cin >> X;
    cout << "Input Alpha: ";
    cin >> Alpha ;
    cout << "Gamma(Alpha,X): " << IncGammaCDF(Alpha, X) << endl;;
  }
	/*  
  else if (Choice==2){
    cout << "Input intial interval X: ";
    cin >> X;
    cout << "Input n: ";
    cin >> n;
    cout << "Input pval: " ;
    cin >> pval;
    cout << "Trunc point" << Zero(X, n/2.0, pval) << endl;
  }
	*/
  else if (Choice==3){
    cout << "X: ";
    cin >> X;
    cout << "m n: ";
    cin >> m >> n;
    cout << "(interval) d1 d2: ";
    cin >> d1 >> d2;
    cout << "j: ";
    cin >> j;
    cout << "g: " << g(X,m,n,d1,d2,Rho,j) << endl;
    cout << "gpp: " << gpp(X, m, n, d1, d2, Rho, j) << endl;

  }
  else if (Choice==4){
    cout << "Input df m,n: " ;
    cin >> m >> n;
    cout << "Input (interval) d1,d2: " ;
    cin >> d1 >> d2;
    cout << "Input (interval) Rho: " ;
    cin >> Rho;
    cout << "Eps1: ";
    cin >> Eps1;
    cout << "Eps2: ";
    cin >> Eps2;
    cout << "nccount: ";
    cin >> nccount;
    cout << "Width: " ;
    cin >> Width;
    BivF(m,n,d1,d2,Rho,Eps1,Eps2,nccount,Width);
  }
  else if (Choice==5){
    cout << "Input interval (a,b): " ;
    cin >> X;
    cout << "Input df m,n: " ;
    cin >> m >> n;
    cout << "Input (interval) d1,d2: " ;
    cin >> d1 >> d2;
    cout << "(interval) Rho: " ;
    cin >> Rho;
    cout << "j: " ;
    cin >> j;
    cout << "nccount: " ;
    cin >> nccount;
    cout << "Newton Cotes: " << NewtonCotes(X,m,n,d1,d2,Rho,j,nccount) << endl;
  }
  else if (Choice==6){
    cout << "Input df m,n: " ;
    cin >> m >> n;
    cout << "Input (interval) d1,d2: " ;
    cin >> d1 >> d2;
    cout << "j: " ;
    cin >> j;
    cout << "(interval) Rho: " ;
    cin >> Rho;
    cout << "Eps1: ";
    cin >> Eps1;
    cout << "Eps2: ";
    cin >> Eps2;
    cout << "nccount: ";
    cin >> nccount;
    cout << "Width: " ;
    cin >> Width;
    cout << "B_j: " << B(m,n,d1,d2,j, Rho,Eps1,Eps2,nccount,Width) << endl;
  }
  else if (Choice==7){
    cout << "x_iminus1, x_i: " << endl;
    cin >> x_iminus1 >> x_i;
    cout << "(interval) Alpha: " ;
    cin >> Alpha;
    cout << "Input df m,n: " ;
    cin >> m >> n;
    cout << "(interval) Rho: " ;
    cin >> Rho;
    cout << "Eps1: ";
    cin >> Eps1;
    cout << "Eps2: ";
    cin >> Eps2;
    cout << "nccount: ";
    cin >> nccount;
    cout << "Width: " ;
    cin >> Width;

    INTERVAL CP, I, S;
    CP=InvBivF(x_iminus1, x_i, Alpha, m, n, Rho, Eps1, Eps2, nccount,Width);
    I = BivF(m,n,Inf(CP),Inf(CP),Rho,Eps1,Eps2,nccount,Width);
    S = BivF(m,n,Sup(CP),Sup(CP),Rho,Eps1,Eps2,nccount,Width);
    cout << endl;
    cout << "CritPt = " <<  CP << endl;
    cout << "Prob   = " << Hull(I,S) << endl;

    cout << endl;
    cout << "x_iminus1: " << x_iminus1 << endl;
    cout << "x_i" << x_i << endl;
    cout << "Alpha: " << Alpha << endl;
    cout << "df m n: " << m << "  " << n << endl; 
    cout << "Rho: " << Rho << endl ;
    cout << "Eps1: " << Eps1 << endl;
    cout << "Eps2: " << Eps2 << endl;
    cout << "nccount: " << nccount << endl;
    cout << "Width: " << Width << endl;

  }
  else if (Choice==8){
    // Find a whole row of critical points in the table
    for (INT i = 2; i+=2; i<4){
      cout << "x_iminus1,x_i: ";
      cin >> x_iminus1 >> x_i;
      cout << "Alpha: ";
      cin >> Alpha;
      cout << "Input df m n: "; 
      cin >> m >> n;
      cout << "(interval) Rho: " ;
      cin >> Rho;
      cout << "Eps1: ";
      cin >> Eps1;
      cout << "Eps2: ";
      cin >> Eps2;
      cout << "nccount: ";
      cin >> nccount;
      cout << "Width: " ;
      cin >> Width;
      
      INTERVAL CP, I, S;
      CP=InvBivF(x_iminus1, x_i, Alpha, m, n, Rho, Eps1, Eps2, nccount,Width);
      I = BivF(m,n,Inf(CP),Inf(CP),Rho,Eps1,Eps2,nccount,Width);
      S = BivF(m,n,Sup(CP),Sup(CP),Rho,Eps1,Eps2,nccount,Width);
      cout << endl;
      cout << "CritPt = " <<  CP << endl;
      cout << "Prob   = " << Hull(I,S) << endl;

      cout << endl;
      cout << "x_iminus1: " << x_iminus1 << endl;
      cout << "x_i" << x_i << endl;
      cout << "Alpha: " << Alpha << endl;
      cout << "Input df m n: " << m << n << endl; 
      cout << "Rho: " << Rho << endl ;
      cout << "Eps1: " << Eps1 << endl;
      cout << "Eps2: " << Eps2 << endl;
      cout << "nccount: " << nccount << endl;
      cout << "Width: " << Width << endl;
    }
  }


}  // end main






