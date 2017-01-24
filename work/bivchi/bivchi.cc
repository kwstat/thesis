// bivchi.cc

// 7.1.98 Changed the p-value to be an INTERVAL
// 6.30.98 Commented out "Ratio" stuff.

#include "MiscFunctions.h"    // For the Random() and Rand01() functions.
#include "Functions.h"      // For Cos, Sqr, Constant::Pi
#include "UnconstrainedOpt.h"
#include "Interval.h"         // for Intersection
#include "Constants.h"        // for Machine::PosInfinity
#include "TestMatrices.h"     // For "matrix" stuff.
#include "Utilities.h"
#include "IntegerVector.h"

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
  // This function 
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
  INTERVAL HAlpha = Hull(Alpha);

  if (Alpha == (int)Alpha){  // Case 1.  Alpha is an integer ----------------

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

    INT n = (int)Alpha;  // n = integer part of alpha
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
    // Split the integral into two pieces.  Use an alternating series for the
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
    INTERVAL EvenSum = Term/( HAlpha+i );       
    
    // OddSum is the sum of the series through an odd term 
    i = 1;
    Term *= -x0/i;       
    INTERVAL OddSum = EvenSum + Term/( HAlpha+i );  // First term i=1
    REAL PrevWidth = 100.0;
    REAL Width = 1.0;
    INT MaxIter = 1000;
    
    // Since this is an alternating series, we calculate two terms at
    // a time and use these as the upper and lower bounds on the sum
    // until the width of the interval doesn't change.
    while ((Width < PrevWidth) && i < MaxIter){             // F
      i++;  // After this line, i is even, so calculate an even term.
      Term *= -x0/i;
      EvenSum = OddSum + Term/( HAlpha+i );
      
      i++;  // After this line, i is odd
      Term *= -x0/i;
      OddSum = EvenSum + Term/( HAlpha+i );
      
      PrevWidth = Width;
      Width = Sup(EvenSum) - Inf(OddSum);
    } 
    INTERVAL AltSeries = Hull(Inf(OddSum), Sup(EvenSum));
    // Now premultiply by:  x0^a / Gamma(a)
    AltSeries *= Exp( HAlpha*Log(x0) - LogGamma(HAlpha) );
    
    // OK, now we compute 1/Gamma(a) * int_x^x0 e^-t t^(a-1) by using
    // a Taylor series.   Left term of RHS of (5.12)
    // First thing we do is pull a trick and change the lower endpoint
    // of the interval 
    x = Hull(Inf(MidPt), Sup(x));
    INTERVAL Taylor = TaylorSeries(x, Md, HAlpha, h);

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

  INTERVAL Result = 0.0
;
  if (Inf(Alpha) == Sup(Alpha))
    Result = IncGammaCDF(Inf(Alpha), X);
  else
    Result = Hull( IncGammaCDF(Inf(Alpha),X),IncGammaCDF(Sup(Alpha),X) );
  return (Result);
}

INTERVAL AdHocJ(INTERVAL c, INTERVAL d, INT m, INT i){
  INTERVAL Result = 0.0;
  INTERVAL Alpha = Hull(m/2.0 + i);

  // The following "if" just makes things run a little faster when c = 0
  if (Sup(c) > 0)
    Result = IncGammaCDF(Alpha, d) - IncGammaCDF(Alpha, c);
  else
    Result = IncGammaCDF(Alpha, d);

  return (Result);
}

INTERVAL Biv1(INTERVAL c1, INTERVAL d1, INTERVAL c2, INTERVAL d2, INT m, INTERVAL Rho){
  // Calculate the density of the bivariate chi-square distribution over
  // a rectangle.  See "Handbook of Statistics", Krishnaiah, p. 754

  INTERVAL J_1i, J_2i, ProbSum = 0.0;
  INTERVAL ErrBoundSum = 0.0, ErrBound = 0.0, PrevErrBound = 0.0, Term;
  INTERVAL Rho2 = Sqr(Rho);
  INTERVAL HalfM = Hull(m/2.0);
  INTERVAL GammaHalfM = Gamma(HalfM);
  INTERVAL ErrBoundConst = Power(1-Rho2,HalfM);

  INTERVAL ScaleFact = 2 * (1-Rho2);
  INTERVAL c1star = c1 / ScaleFact, d1star = d1 / ScaleFact;
  INTERVAL c2star = c2 / ScaleFact, d2star = d2 / ScaleFact;

  BOOL Converge = 0;
  INT i = 0;
  INTERVAL Numer = 1.0, Denom = 1.0;
  while (!Converge){
    // This is just hoping for a little more speed...about 5%
    if ((c1star==c2star) && (d1star==d2star))
      J_1i = J_2i = AdHocJ(c1star, d1star, m, i);
    else{
      J_1i = AdHocJ(c1star, d1star, m, i);
      J_2i = AdHocJ(c2star, d2star, m, i);
    }
    // I made substantial improvements in accuracy by eliminating the calls
    // to the Gamma function
    Term = Numer / Denom;
    ProbSum += Term * J_1i * J_2i;
    ErrBoundSum += Term;
    ErrBound = 1-ErrBoundConst * ErrBoundSum;

    if (PrevErrBound == ErrBound) 
      Converge = 1;
    else{
      // Now prepare for next iteration of while loop
      i++;
      // Rescale from time to time to prevent overflow
      if (Sup(Numer) > 1.0e300 || Sup(Denom) > 1.0e300) {
        Numer /= 1.0e150;
        Denom /= 1.0e150;
      }
      Numer *= (m + 2*(i-1)) * Rho2;
      Denom *= i * 2;
      PrevErrBound = ErrBound;
    }
  }

  // First adjust for the constant multiplier to get the probability
  INTERVAL Prob = ProbSum * Power(1-Rho2,HalfM);
  // Then enlarge the right endpoint by the error bound
  Prob = Hull(Inf(Prob), Sup(Prob)+ErrBound);

  return (Prob);
}

INTERVAL InvBiv1(REAL x_iminus1, REAL x_i, INTERVAL Prob, 
                       INT m, INTERVAL Rho){
  INT i=0;  
  BOOL Done = 0;
  REAL L = x_iminus1;
  REAL U = x_i;
  INTERVAL FL = Biv1(0, L, 0, L, m, Rho) - Prob;
  INTERVAL FU = Biv1(0, U, 0, U, m, Rho) - Prob;

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
    F_iplus1 = Biv1(0, x_iplus1, 0, x_iplus1, m, Rho) - Prob;

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
  cout << "Bracket secant : " << i << " " << Hull(L,U) << endl;
 
  // Now switch to bisection search
  // First narrow the interval by moving the right endpoint
  PrevQ = PrevU = -1;
  Q = L;
  while ( (Q != PrevQ) || (U != PrevU) ) {
    PrevQ = Q;
    PrevU = U;
    i++;
    x = (Q + U) / 2.0;
    Fx = Biv1(0, x, 0, x, m, Rho) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) Q = L = x;
    else Q = x;
  }
  cout << "Right bisection: " << i << " " << Hull(L,U) << endl;

  // Then narrow the interval by moving the left endpoint
  PrevQ = PrevL = -1;
  Q = U;
  while ( (Q != PrevQ) || (L != PrevL) ) {
    PrevQ = Q;
    PrevL = L;
    i++;
    x = (L + Q) / 2.0;
    Fx = Biv1(0, x, 0, x, m, Rho) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) L = x;
    else Q = x;
  }
  cout << "Left  bisection: " << i << " " << Hull(L,U) << endl;

  return (Hull(L,U));
}

INTERVAL Biv2(INTERVAL c1, INTERVAL c2, INT m, INT n, INTERVAL Rho){
  INTERVAL J_i, K_ij, ProbSum = 0.0;
  INTERVAL ErrBoundSum = 0.0, ErrBound = 0.0;
  INTERVAL iPrevErrBoundSum, jPrevErrBoundSum;
  INTERVAL Term;
  INTERVAL Rho2 = Sqr(Rho);
  INTERVAL HalfM = Hull(m/2.0), HalfN = Hull(n/2.0);
  INTERVAL GammaHalfM = Gamma(HalfM), GammaHalfN = Gamma(HalfN);
  INTERVAL ErrBoundConst = Power(1-Rho2,HalfM+HalfN);

  INTERVAL ScaleFact = 2 * (1-Rho2);
  INTERVAL c1star = c1 / ScaleFact;
  INTERVAL c2star = c2 / ScaleFact;

  BOOL iConverge = 0, jConverge = 0;
  INT i = 0, j = 0;
  INTERVAL iNumer, iDenom, iTerm;
  INTERVAL jNumer, jDenom, jTerm;

  iConverge = 0;
  i = 0;
  iPrevErrBoundSum = 0.0;
  iNumer = iDenom = 1.0;
  while (!iConverge){
    J_i = IncGammaCDF(Hull(m/2 +i), c1star);
    iTerm = iNumer / iDenom;
  
    jConverge = 0;
    j = 0;
    jPrevErrBoundSum = 0.0;
    jNumer = jDenom = 1.0;
    // ----------------------------------------------------------------
    while (!jConverge){
      K_ij = IncGammaCDF(Hull((m+n)/2 +i+j), c2star);
      jTerm = jNumer / jDenom;
      Term = iTerm * jTerm;
      ProbSum += Term * J_i * K_ij;
      ErrBoundSum += Term;
      ErrBound = 1 - ErrBoundConst * ErrBoundSum;

      j++;
      if (Inf(ErrBoundSum) == Inf(jPrevErrBoundSum) ) jConverge = 1;
      else {   // Prepare for next iteration of j loop
        // Rescale from time to time to prevent overflow
        if (Sup(jNumer) > 1.0e300 || Sup(jDenom) > 1.0e300) {
          jNumer /= 1.0e150;
          jDenom /= 1.0e150;
        }
        jNumer *= (n + 2*(j-1)) * Rho2;
        jDenom *= j * 2;
        jPrevErrBoundSum = ErrBoundSum;
      }

    } //-----------------------------------------------------------------

    i++;
    if (Inf(ErrBoundSum) == Inf(iPrevErrBoundSum) ) iConverge = 1;
    else {  // Prepare for next iteration of i loop
      // Rescale from time to time to prevent overflow
      if (Sup(iNumer) > 1.0e300 || Sup(iDenom) > 1.0e300) {
        iNumer /= 1.0e150;
        iDenom /= 1.0e150;
      }
      iNumer *= (m + 2*(i-1)) * Rho2;
      iDenom *= i * 2;
      iPrevErrBoundSum = ErrBoundSum;
    }

  }
  // First adjust for the constant multiplier to get the probability
  INTERVAL Prob = ProbSum * ErrBoundConst;
  // Then enlarge the right endpoint by the error bound
  Prob = Hull(Prob, Prob+ErrBound);
  return (Prob);
}

INTERVAL InvBiv2(REAL x_iminus1, REAL x_i, INTERVAL Prob, 
                 INT m, INT n, INTERVAL Rho){
  //REAL Ratio;
  //cout << "Input interval ratio C1/C2:" ;
  //cin >> Ratio;
  //cout << endl;

  INT i=0;  
  BOOL Done = 0;
  REAL L = x_iminus1;
  REAL U = x_i;
  //INTERVAL FL = Biv2(Ratio*L, L, m, n, Rho) - Prob;
  //INTERVAL FU = Biv2(Ratio*U, U, m, n, Rho) - Prob;
  INTERVAL FL = Biv2(L, L, m, n, Rho) - Prob;
  INTERVAL FU = Biv2(U, U, m, n, Rho) - Prob;
  INTERVAL F_i = Biv2(x_i, x_i, m, n, Rho) - Prob;

  cout << "FL " << FL << endl;
  cout << "FU " << FU << endl;
 
  REAL Q, PrevQ, PrevL, PrevU, x, x_iplus1;
  INTERVAL Fx, F_iplus1, IX_c, IX_i, IX_iplus1, X_iplus1;

  if (Sup(FL) < 0  &&  Inf(FU) > 0) ;
    // values are satisfactory
  else {
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
    //F_iplus1 = Biv2(Ratio*x_iplus1, x_iplus1, m, n, Rho) - Prob;
    F_iplus1 = Biv2(x_iplus1, x_iplus1, m, n, Rho) - Prob;

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
  cout << "Bracket secant : " << i << " " << Hull(L,U) << endl;
 
  // Now switch to bisection search
  // First narrow the interval by moving the right endpoint
  PrevQ = PrevU = -1;
  Q = L;
  while ( (Q != PrevQ) || (U != PrevU) ) {
    PrevQ = Q;
    PrevU = U;
    i++;
    x = (Q + U) / 2.0;
    //Fx = Biv2(Ratio*x, x, m, n, Rho) - Prob;
    Fx = Biv2(x, x, m, n, Rho) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) Q = L = x;
    else Q = x;
  }
  cout << "Right bisection: " << i << " " << Hull(L,U) << endl;

  // Then narrow the interval by moving the left endpoint
  PrevQ = PrevL = -1;
  Q = U;
  while ( (Q != PrevQ) || (L != PrevL) ) {
    PrevQ = Q;
    PrevL = L;
    i++;
    x = (L + Q) / 2.0;
    //Fx = Biv2(Ratio*x, x, m, n, Rho) - Prob;
    Fx = Biv2(x, x, m, n, Rho) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) L = x;
    else Q = x;
  }
  cout << "Left  bisection: " << i << " " << Hull(L,U) << endl;

  //if (Ratio != 1.0){
  //  cout << "C1 = " << Ratio*Hull(L,U) << endl;
  //  cout << "C2 = " << Hull(L,U) << endl;
  //}

  return (Hull(L,U));
}

INTERVAL Biv3(INTERVAL c1, INTERVAL c2, INT m, INT n, INT p, INTERVAL Rho){

  INTERVAL J_jk, K_jl, ProbSum = 0.0;
  INTERVAL ErrBoundSum = 0.0, ErrBound = 0.0;
  INTERVAL jPrevErrBoundSum, kPrevErrBoundSum, lPrevErrBoundSum;
  INTERVAL Term;
  INTERVAL Rho2 = Sqr(Rho);
  INTERVAL HalfM = Hull(m/2.0), GammaHalfM = Gamma(HalfM);
  INTERVAL HalfN = Hull(n/2.0), GammaHalfN = Gamma(HalfN);
  INTERVAL HalfP = Hull(p/2.0), GammaHalfP = Gamma(HalfP);
  INTERVAL ErrBoundConst = Power(1-Rho2,HalfM+HalfN+HalfP);

  INTERVAL ScaleFact = 2 * (1-Rho2);
  INTERVAL c1star = c1 / ScaleFact;
  INTERVAL c2star = c2 / ScaleFact;

  BOOL jConverge = 0, kConverge = 0, lConverge = 0;
  INT j = 0, k = 0, l = 0;
  INTERVAL jNumer, jDenom, jTerm;
  INTERVAL kNumer, kDenom, kTerm;
  INTERVAL lNumer, lDenom, lTerm;

  jConverge = 0;
  j = 0;
  jPrevErrBoundSum = 0.0;
  jNumer = jDenom = 1.0;
  while (!jConverge){  //---------------------------------------------------
    jTerm = jNumer / jDenom;
  
    kConverge = 0;
    k = 0;
    kPrevErrBoundSum = 0.0;
    kNumer = kDenom = 1.0;
    while (!kConverge){  //-------------------------------------------------
      J_jk = IncGammaCDF(Hull((m+n)/2 +j +k), c1star);
      kTerm = kNumer / kDenom;

      lConverge = 0;
      l = 0;
      lPrevErrBoundSum = 0.0;
      lNumer = lDenom = 1.0;
      while (!lConverge){  //-----------------------------------------------
        K_jl = IncGammaCDF(Hull((m+p)/2 +j +l), c2star);
        lTerm = lNumer / lDenom;
        Term = jTerm * kTerm * lTerm;
        ProbSum += Term * J_jk * K_jl;
        ErrBoundSum += Term;
        ErrBound = 1 - ErrBoundConst * ErrBoundSum;

        l++;
        if (Inf(ErrBoundSum) == Inf(lPrevErrBoundSum) ) lConverge = 1;
        else { // Prepare for next iteration of l loop
          if (Sup(lNumer) > 1.0e300 || Sup(lDenom) > 1.0e300) {
            lNumer /= 1.0e150;
            lDenom /= 1.0e150;
          }
          lNumer *= (p + 2*(l-1)) * Rho2;
          lDenom *= l * 2;
          lPrevErrBoundSum = ErrBoundSum;
        }
      }  // end l while

      k++;
      if (Inf(ErrBoundSum) == Inf(kPrevErrBoundSum) ) kConverge = 1;
      else {   // Prepare for next iteration of k loop
        // Rescale from time to time to prevent overflow
        if (Sup(kNumer) > 1.0e300 || Sup(kDenom) > 1.0e300) {
          kNumer /= 1.0e150;
          kDenom /= 1.0e150;
        }
        kNumer *= (n + 2*(k-1)) * Rho2;
        kDenom *= k * 2;
        kPrevErrBoundSum = ErrBoundSum;
      }
    } // end k while

    j++;
    if (Inf(ErrBoundSum) == Inf(jPrevErrBoundSum) ) jConverge = 1;
    else {  // Prepare for next iteration of j loop
      // Rescale from time to time to prevent overflow
      if (Sup(jNumer) > 1.0e300 || Sup(jDenom) > 1.0e300) {
        jNumer /= 1.0e150;
        jDenom /= 1.0e150;
      }
      jNumer *= (m + 2*(j-1)) * Rho2;
      jDenom *= j * 2;
      jPrevErrBoundSum = ErrBoundSum;
    }

  }
  // First adjust for the constant multiplier to get the probability
  INTERVAL Prob = ProbSum * ErrBoundConst;
  // Then enlarge the right endpoint by the error bound
  Prob = Hull(Prob, Prob+ErrBound);
  return (Prob);

}

INTERVAL InvBiv3(REAL x_iminus1, REAL x_i, INTERVAL Prob, 
                 INT m, INT n, INT p, INTERVAL Rho){
  INT i=0;  
  BOOL Done = 0;
  REAL L = x_iminus1;
  REAL U = x_i;
  INTERVAL FL = Biv3(L, L, m, n, p, Rho) - Prob;
  INTERVAL FU = Biv3(U, U, m, n, p, Rho) - Prob;
  INTERVAL F_i = Biv3(x_i, x_i, m, n, p, Rho) - Prob;

  cout << "FL " << FL << endl;
  cout << "FU " << FU << endl;
 
  REAL Q, PrevQ, PrevL, PrevU, x, x_iplus1;
  INTERVAL Fx, F_iplus1, IX_c, IX_i, IX_iplus1, X_iplus1;

  if (Sup(FL) < 0  &&  Inf(FU) > 0) ;
    // values are satisfactory
  else {
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
    F_iplus1 = Biv3(x_iplus1, x_iplus1, m, n, p, Rho) - Prob;

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
  cout << "Bracket secant : " << i << " " << Hull(L,U) << endl;
 
  // Now switch to bisection search
  // First narrow the interval by moving the right endpoint
  PrevQ = PrevU = -1;
  Q = L;
  while ( (Q != PrevQ) || (U != PrevU) ) {
    PrevQ = Q;
    PrevU = U;
    i++;
    x = (Q + U) / 2.0;
    Fx = Biv3(x, x, m, n, p, Rho) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) Q = L = x;
    else Q = x;
  }
  cout << "Right bisection: " << i << " " << Hull(L,U) << endl;

  // Then narrow the interval by moving the left endpoint
  PrevQ = PrevL = -1;
  Q = U;
  while ( (Q != PrevQ) || (L != PrevL) ) {
    PrevQ = Q;
    PrevL = L;
    i++;
    x = (L + Q) / 2.0;
    Fx = Biv3(x, x, m, n, p, Rho) - Prob;
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) L = x;
    else Q = x;
  }
  cout << "Left  bisection: " << i << " " << Hull(L,U) << endl;

  return (Hull(L,U));
}


INT main(){
  cout.precision(15);
  INT Choice;

  cout << endl;
  cout << "Bivariate Chi-square" << endl << endl;

  cout << "1  Evaluate Inc Gamma: " << endl;
  cout << endl;
  cout << "Case I: equal df" << endl;
  cout << "2  Biv1     (Cumulative prob.) " << endl;
  cout << "3  InvBiv1  (Critical points ) " << endl;
  cout << "4  InvBiv1  (Multiple rho) " << endl;
  cout << endl;

  cout << "Case II: Unequal df" << endl;
  cout << "5  Biv2 " << endl;
  cout << "6  InvBiv2 " << endl;
  cout << endl;

  cout << "Case III: Unequal df, nonzero correlations" << endl;
  cout << "9  Biv3 " << endl;
  cout << "10 InvBiv3 " << endl;

  cout << "Enter choice: " ;
  cin >> Choice;
  cout << endl;

  INTERVAL Alpha, X, Rho, c1, c2, d1, d2;
  INTERVAL Prob, CritPt, I, S;
  INT m, n, p;
  REAL x0, x1;

  if (Choice==1){
    cout << "Input (interval) Alpha: " ;
    cin >> Alpha;
    cout << "Input X: " ;
    cin >> X;
    cout << "IncGammaCDF(Alpha,X): " << IncGammaCDF(Alpha, X) << endl;
  }
  if (Choice==2){
    cout << "(degenerate interval) c1 = " ;
    cin >> c1;
    cout << "(degenerate interval) d1 = " ;
    cin >> d1;
    cout << "(degenerate interval) c2 = " ;
    cin >> c2;
    cout << "(degenerate interval) d2 = " ;
    cin >> d2;
    cout << "m  = " ;
    cin >> m;
    cout << "(interval) rho = ";
    cin >> Rho;

    Prob = Biv1(c1, d1, c2, d2, m, Rho);
    cout << "Prob = " << Prob << endl;
  }
  
  else if (Choice == 3 || Choice == 4){
    cout << "Input (interval) Prob = P(y1 < p, y2 < p): ";
    cin >> Prob;
    cout << "Input m " ;
    cin >> m;
    if (Choice == 3){
    cout << "Input (interval) Rho " ;
    cin >> Rho;
    }
    cout << "x0 " ;
    cin >> x0;
    cout << "x1 " ;
    cin >> x1;
    cout << "Critical Values : " ;
    cout << "m = " << m << endl;
    if (Choice ==  3){
        CritPt = InvBiv1(x0, x1, Prob, m, Rho);
        // Biv1 expects REAL points, so we evaluate it at the upper and lower
        // ends of CritPt and then take the hull of the two results
        I = Biv1(0, Inf(CritPt), 0, Inf(CritPt), m, Rho);
        S = Biv1(0, Sup(CritPt), 0, Sup(CritPt), m, Rho);
        cout << endl;
        cout << "Rho    = " << Rho << endl;
        cout << "CritPt = " <<  CritPt << endl;
        cout << "Prob   = " << Hull(I,S) << endl;
    }
    else {
      INT kw;
      for (kw = 1; kw < 10; kw += 1){
				Rho = Hull(kw/10.0);
        CritPt = InvBiv1(x0, x1, Prob, m, Rho);
        // Biv1 expects REAL points, so we evaluate it at the upper and lower
        // ends of CritPt and then take the hull of the two results
        I = Biv1(0, Inf(CritPt), 0, Inf(CritPt), m, Rho);
        S = Biv1(0, Sup(CritPt), 0, Sup(CritPt), m, Rho);
        cout << endl;
        cout << "Rho    = " << Rho << endl;
        cout << "CritPt = " <<  CritPt << endl;
        cout << "Prob   = " << Hull(I,S) << endl;
      } // end for
    } // end else choice = 4
    
  }  // end choice 3

  else if (Choice == 5 || Choice == 9) {
    cout << "(degenerate interval) d1 = " ;
    cin >> d1;
    cout << "(degenerate interval) d2 = " ;
    cin >> d2;
    cout << "m  = " ;
    cin >> m;
    cout << "n  = " ;
    cin >> n;
    if (Choice == 9){
      cout << "p = " ;
      cin >> p;
    }
    cout << "(interval) rho = ";
    cin >> Rho;

    if (Choice ==5) Prob = Biv2(d1, d2, m, n, Rho);
    else Prob = Biv3(d1, d2, m, n, p, Rho);

    cout << "Prob = " << Prob << endl;

  }
/*-----------------------------------------------------------------*/


  else if (Choice == 6){
    cout << "Input (interval) Prob P(y1 < p, y2 < p): ";
    cin >> Prob;
    cout << "Input deg freedom m n: " ;
    cin >> m >> n;
    cout << "Input correlation (interval) rho: ";
    cin >> Rho;
    cout << "Input initial brackets x0 x1:" ;
    cin >> x0 >> x1;

    CritPt = InvBiv2(x0, x1, Prob, m, n, Rho);
    // Biv3 expects REAL points, so we evaluate it at the upper and lower
    // ends of CritPt and then take the hull of the two results
    I = Biv2(Inf(CritPt), Inf(CritPt), m, n, Rho);
    S = Biv2(Sup(CritPt), Sup(CritPt), m, n, Rho);
    cout << "CritPt = " <<  CritPt << endl;
    cout << "Prob   = " << Hull(I,S) << endl;
  }  // end choice 6

  else if (Choice == 10){
    cout << "Input (interval) Prob = P(y1 < p, y2 < p): ";
    cin >> Prob;
    cout << "Input m n p:" ;
    cin >> m >> n >> p;
    cout << "Input (interval) Rho: ";
    cin >> Rho;
    cout << "Input x0 x1:" ;
    cin >> x0 >> x1;

    CritPt = InvBiv3(x0, x1, Prob, m, n, p, Rho);
    // Biv3 expects REAL points, so we evaluate it at the upper and lower
    // ends of CritPt and then take the hull of the two results
    I = Biv3(Inf(CritPt), Inf(CritPt), m, n, p, Rho);
    S = Biv3(Sup(CritPt), Sup(CritPt), m, n, p, Rho);
    cout << "CritPt = " <<  CritPt << endl;
    cout << "Prob   = " << Hull(I,S) << endl;
  }  // end choice 10

}  // end main
