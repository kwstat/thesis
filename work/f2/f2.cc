// f2.cc

#include "MiscFunctions.h"    // For the Random() and Rand01() functions.
#include "Functions.h"      // For Cos, Sqr, Constant::Pi
#include "UnconstrainedOpt.h"
#include "Interval.h"         // for Intersection
#include "Constants.h"        // for Machine::PosInfinity
#include "Utilities.h"
#include "IntegerVector.h"

INTERVAL Comb(INT m, INT n){
  // Compute  m!/[n!(m-n)!]

  INTERVAL Result = 1.0;
  INT k;

  if (n < m-n) 
    n = m-n;

  for(k=m;k>n; k--){
		Result *= k;
	  Result /= k-n;
	}
  return (Result);
}

INTERVAL q(INTERVAL H, INT n, INT m){
  // q(H), See Hamdy, page 160.
  INT i, j;
  INTERVAL Result = 0.0;

  for(i=0; i<n; i++)
		for(j=0; j<m; j++){
			Result += Comb(n+i-1,i)*Comb(n+m+i-1,j)*Power(H,j)*
                  Power(1-H,n+m+i-1-j)/Power(2.0,n+i-1);
		}
  return (Result);
}

INTERVAL Invq(REAL x_iminus1, REAL x_i, REAL pval, INT m, INT n){
  // Note, the code was adapted from bivchi.cc, which assumed that the
  // CDF is increasing.  Since the q() function is decreasing, I threw in
  // -q() everywhere.

  INT i=0;  
  BOOL Done = 0;
  REAL L = x_iminus1;
  REAL U = x_i;
  INTERVAL FL = -(1-q(Hull(x_iminus1), m, n) - pval);
  INTERVAL FU = -(1-q(Hull(x_i), m, n) - pval);

  INTERVAL F_i = -(1-q(Hull(x_i), m, n) - pval);
 
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
		F_iplus1 = -(1-q(Hull(x_iplus1), m, n) - pval);

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

  // Now switch to bisection search
  // First narrow the interval by moving the right endpoint
  PrevQ = PrevU = -1;
  Q = L;
  while ( (Q != PrevQ) || (U != PrevU) ) {
    PrevQ = Q;
    PrevU = U;
    i++;
    x = (Q + U) / 2.0;
    Fx = -(1-q(Hull(x), m, n) - pval);
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) Q = L = x;
    else Q = x;
  }

  // Then narrow the interval by moving the left endpoint
  PrevQ = PrevL = -1;
  Q = U;
  while ( (Q != PrevQ) || (L != PrevL) ) {
    PrevQ = Q;
    PrevL = L;
    i++;
    x = (L + Q) / 2.0;
    Fx = -(1-q(Hull(x), m, n) - pval);
    if (Inf(Fx) > 0) U = x;
    else if (Sup(Fx) < 0) L = x;
    else Q = x;
  }

  return (Hull(L,U));
}

INT main(){
  cout.precision(16);
  INT Choice;

  REAL x0, x1, pval;
  INT m,n;
  REAL h, c;
  INTERVAL H, C;

  cout << "1. CDF" << endl;
  cout << "2. CDF points" << endl;
  cout << "3. InvCDF" << endl;
  cout << "4. Reproduce Gupta/Sobel table for df n := m := nu" << endl;
  cout << "Input choice: " ;
  cin >> Choice;

  if (Choice==1){
		cout << "Input INT deg. freedom m n: " ;
		cin >> m >> n;    
    cout << "Input c: ";
    cin >> c;
    h = 1/(1+2*n*c/m);
    cout << "Prob: " << 1-q(Hull(h), m, n) << endl;;    
	}
  else if (Choice ==2){
    cout << "Input INT deg. freedom m n: " ;
		cin >> m >> n;    
    for (c=0; c<=1; c+= .1){
			h = 1/(1+2*n*c/m);
			cout <<"c: " << c <<" h: " <<h <<" Prob: " <<1-q(Hull(h), m, n) << endl;
		}    
  }
  else if (Choice ==3){
		cout << "Input REAL P(Y < c) = pval: " ;
		cin >> pval;
		cout << "Input INT deg. freedom m n: " ;
		cin >> m >> n;
		cout << "Input bracketing REAL x0 x1: ";
		cin >> x0 >> x1;
    // translate from c to h
    x0 = 1/(1+2*n*x0/m);
    x1 = 1/(1+2*n*x1/m);
      cout << "x0 " << x0 << endl;
			cout << "x1 " << x1 << endl;
			cout << "pval: " << pval << endl;
			cout << "m: " << m << endl;
			cout << "n: " << n << endl;

    H = Invq(x0, x1, pval, m, n);
		// then translate back to c
    C = m*(1/H -1) / (2.0*n);
		cout << "Crit Pt: " << C << endl;
  }
  else if (Choice ==4){
    pval = 0.75;
		cout << endl;
		cout << "pval: " << pval << endl;
		for (m=1; m<= 25; m++){
			n = m;
			x0 = .001;
			x1 = 1.0;
			x0 = 1/(1+2*n*x0/m);
			x1 = 1/(1+2*n*x1/m);
			H = Invq(x0, x1, pval, m, n);
			// then translate back to c
			C = m*(1/H -1) / (2.0*n);
			cout << "df m = n:" << m << " Crit Pt: " << C << endl;
		}

    pval = 0.90;
		cout << endl;
		cout << "pval: " << pval << endl;
		for (m=1; m<= 25; m++){
			n = m;
			x0 = .001;
			x1 = 1.0;
			x0 = 1/(1+2*n*x0/m);
			x1 = 1/(1+2*n*x1/m);
			H = Invq(x0, x1, pval, m, n);
			// then translate back to c
			C = m*(1/H -1) / (2.0*n);
			cout << "df m = n:" << m << " Crit Pt: " << C << endl;
		}

    pval = 0.95;
		cout << endl;
		cout << "pval: " << pval << endl;
		for (m=1; m<= 25; m++){
			n = m;
			x0 = .001;
			x1 = 1.0;
			x0 = 1/(1+2*n*x0/m);
			x1 = 1/(1+2*n*x1/m);
			H = Invq(x0, x1, pval, m, n);
			// then translate back to c
			C = m*(1/H -1) / (2.0*n);
			cout << "df m = n:" << m << " Crit Pt: " << C << endl;
		}

    pval = 0.99;
		cout << endl;
		cout << "pval: " << pval << endl;
		for (m=1; m<= 25; m++){
			n = m;
			x0 = .001;
			x1 = 1.0;
			x0 = 1/(1+2*n*x0/m);
			x1 = 1/(1+2*n*x1/m);
			H = Invq(x0, x1, pval, m, n);
			// then translate back to c
			C = m*(1/H -1) / (2.0*n);
			cout << "df m = n:" << m << " Crit Pt: " << C << endl;
		}
	} // End of choice 4
}


