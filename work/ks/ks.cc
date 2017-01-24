// ks.cc
// Distribution function of the Kolmogorov-Smirnov statistic

#include "Interval.h"
#include "Functions.h"  // Exp

INTERVAL KSDens(REAL x){
  INTERVAL OddSum = 0, EvenSum = 0;
	INTERVAL Sum = 0;
  INTERVAL CumEnc = Hull(-10,10), PrevEnc = Hull(-10,10);

  BOOL Converge = 0;
  INT j = 0;

  while (!Converge){
    PrevEnc = CumEnc;

		j++;
    OddSum = EvenSum + Sqr(j) * Exp(-2*Sqr(j*x));
		j++;
		EvenSum = OddSum - Sqr(j) * Exp(-2*Sqr(j*x));
	  Sum = Hull(OddSum, EvenSum);

    Intersection(CumEnc, PrevEnc, Sum);

//		cout << "j: " << j << "   CumEnc: " << CumEnc << endl;
		if (CumEnc == PrevEnc) Converge = 1;
	}
  return (8 * x * Sum);

}

INTERVAL KSDist(REAL x){
  INTERVAL OddSum = 0, EvenSum = 0;
	INTERVAL Sum = 0;
  INTERVAL CumEnc = Hull(-10,10), PrevEnc = Hull(-10,10);

  BOOL Converge = 0;
  INT j = 0;

  while (!Converge){
    PrevEnc = CumEnc;

		j++;
    OddSum = EvenSum + Exp(-2*Sqr(j*x));
		j++;
		EvenSum = OddSum - Exp(-2*Sqr(j*x));
	  Sum = Hull(OddSum, EvenSum);

    Intersection(CumEnc, PrevEnc, Sum);

//		cout << "j: " << j << "   CumEnc: " << CumEnc << endl;
		if (CumEnc == PrevEnc) Converge = 1;
	}
  return (1 - 2* Sum);
}

INTERVAL f(INTERVAL x, REAL alpha){
  INTERVAL R;
  R = Hull(KSDist(Inf(x))-alpha, KSDist(Sup(x))-alpha);
  return (R);
}

INTERVAL f(REAL x, REAL alpha){
  INTERVAL R;
  R = KSDist(x) - alpha;
  return (R);
}

INTERVAL fp(INTERVAL x){
  INTERVAL R;
  R = Hull(KSDens(Inf(x)), KSDens(Sup(x)));
  return (R);
}

INTERVAL Zero(REAL alpha, INTERVAL x){
  // This function needs f(x,alpha) and fp(X)

  BOOL Done = 0;
  INTERVAL Xn = x, Xnp1, Xtemp;
  REAL xm;

  if (0 <= fp(Xn)) Done = 1;
	if (Sup( f(Inf(Xn),alpha)*f(Sup(Xn),alpha) )>= 0) Done = 1;

  while (!Done){
    
    // Now compute the next iterate
    xm = (Inf(Xn)+Sup(Xn)) / 2.0;
    Xtemp = xm - f(xm, alpha) / fp(Xn);

		Intersection(Xnp1, Xn, Xtemp);

		if (Xn == Xnp1) Done = 1;

    // Prepare for next iteration									 
    Xn = Xnp1;
		if (0 <= fp(Xn)) Done = 1;
		if (Sup( f(Inf(Xn),alpha)*f(Sup(Xn),alpha) )>= 0) Done = 1;
	}

  return (Xn);
}

VOID main(){
  REAL alpha, x;
  INTERVAL X0;
  INT Choice;

  cout << "1. F(x) " << endl;
  cout << "2. Crit Pt " << endl;
  cout << "3. Check tables (15 digits)" << endl;

  cout <<  "Choice: " << endl;
  cin >> Choice;

  if (Choice==1) {
    cout.precision(16);
    cout << "Enclosure of Kolmogorov-Smirnov distribution for real x." << endl;
		cout << "Input x: " ;
		cin >> x;
		cout << endl;
		cout << "f(x) = " << KSDens(x) << endl;
		cout << "F(x) = " << KSDist(x) << endl;
	}
  else if (Choice==2){
		cout << "Interval-Newton iteration for finding percentile point." << endl; 
		cout << "Percentile: " ;
		cin >> alpha;
		cout << endl;
		cout << "Initial zero-bracketing interval X0: " ;
		cin >> X0;
		cout << endl;
		cout.precision(17);
		cout << "Zero(alpha): " << Zero(alpha, X0) << endl;
  }
  else if (Choice==3){
		cout << "Check of tables in Smirnov (1948)" << endl;
		cout.precision(16);
		cout << endl;
		for(x=0.28; x <= 3.00; x+=.01){
			cout << "x: " << x << "  F(x): " << KSDist(x) << endl;
		}
	}

}
