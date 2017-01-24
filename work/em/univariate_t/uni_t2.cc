// uni_t2.cc
// Psi = Mu
// Univariate t example with multiple stationary points.
// Example from EM book

#include "UnconstrainedOpt.h"
#include "Functions.h"
//#define DEBUG


// Globals
INTERVAL Psi_k;
INTERVAL_VECTOR W(6);
INTERVAL Nu;



INTERVAL GradientQ(INTERVAL_VECTOR & Psi)
{
  INTERVAL Mu = Psi(1);
  INTERVAL Result = 0;
  INT i, lenw=Dimension(W);

  for(i=1;i<=lenw;i++){
		Result += 2*(W(i)-Mu)/(Nu+ Sqr(W(i)-Psi_k) ) ;
	}

  return (Result);
}

VOID SetTestDomain (INTERVAL_VECTOR &Psi)
{
  Resize(Psi, 1);
  cout << "Enter interval to enclose scalar Psi (Mu): ";
  cin >> Psi(1);
}

INTERVAL IntervalFun (INTERVAL &Psi){
  INTERVAL Mu = Psi;
  INT i, lenw = Dimension(W);
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
    F += Sqr(W(i)-Mu) / (Nu+Sqr(W(i)-Psi_k) );
  F *= -(Nu + 1);
  return(F);
}

REAL RealUpper (VECTOR & Psi)
// real valued test function
{
  REAL Mu = Psi(1);
  INT i, lenw = Dimension(W);
  REAL f = 0;
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
    F += Sqr(W(i)-Mu) / (Nu+Sqr(W(i)-Psi_k) );
  F *= (Nu + 1);
  f = Sup(F);

#ifdef DEBUG
  cout << "Real f(Theta):  " << f << endl;
#endif 
  return (f);
}

INTERVAL IntervalUpper (INTERVAL_VECTOR & Psi)
// interval valued test function
{
  INTERVAL Mu = Psi(1);
  INT i, lenw = Dimension(W);
  REAL f = 0;
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
		F += Sqr((W(i)-Mu))/(Nu + Sqr((W(i)-Psi_k)));
  F *= (Nu + 1);
  f = Sup(F);

  F = Hull(f);

#ifdef DEBUG
  cout << "Interval f(Theta): " << F << endl;
#endif
  return (F);
}


REAL RealLower (VECTOR &Psi)
// real valued test function
{
  REAL Mu = Psi(1);
  INT i, lenw = Dimension(W);
  REAL f = 0;
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
		F += Sqr((W(i)-Mu))/(Nu + Sqr((W(i)-Psi_k)));
  F *= -(Nu + 1);
  f = Inf(F);

#ifdef DEBUG
  cout << "Real f(Theta):  " << f << endl;
#endif 
  return (f);
}

INTERVAL IntervalLower (INTERVAL_VECTOR &Psi)
// interval valued test function
{
  INTERVAL Mu = Psi(1);
  INT i, lenw = Dimension(W);
  REAL f = 0;
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
		F += Sqr((W(i)-Mu))/(Nu + Sqr((W(i)-Psi_k)));
  F *= -(Nu + 1);
  f = Inf(F);

  F = Hull(f);

#ifdef DEBUG
  cout << "Interval f(Theta): " << F << endl;
#endif
  return (F);
}


INT main()
{
  INT LenW;
  // First define the data
  cout << "Dimension of W: " ;
  cin >> LenW;
  Resize(W,LenW);
  cout << "Enter real vector W: ";
  cin >> W;
  //  W(1) = -20;  W(2) = 1; W(3) = 2; W(4) = 3;

  // cout << "Enter interval Nu: ";
  // cin >> Nu;
  Nu = Hull(.05);

  INT Iterations;
  INTERVAL_VECTOR TestDomain, PointDomain(1), PointU, PointL;
  INTERVAL_VECTOR IVPsi_k;
  REAL LowerBound, UpperBound;

  INT Narrow, i, EMiter, Method, IntersectYN;
  INTERVAL PrevPsi_k, Psi_kL, Psi_kR;

  cout << "Enter number of EM iterations EMiter: " ;
  cin >> EMiter;

  cout << "Enter interval Psi_k (i.e. M_0^t)";
  cin >> Psi_k;
  INTERVAL OrigPsi = Psi_k;

	SetTestDomain(TestDomain);


  INT PreSrch;
  REAL BoxWdth, BoxLeft;
  INTERVAL_VECTOR PreBox(1);
  cout << "Enter 0 to skip the presearch" << endl;
  cout << "Enter 1 to do the pre-search." << endl;
  cin >> PreSrch;
  if (PreSrch==1){

    cout << "Enter pre-search box width" << endl;
    cin >> BoxWdth;

    REAL SupPsi = Sup(Psi_k);
    for (BoxLeft=Inf(Psi_k); BoxLeft<SupPsi; 
           BoxLeft = Inf(Hull(BoxLeft+BoxWdth))) {
      PreBox(1) = Hull(BoxLeft,Sup(Hull(BoxLeft+BoxWdth)));
      Psi_k=PreBox(1);
      if (!(0 <= GradientQ(PreBox)))
        cout << "No   . " ;
			else cout << "Maybe. " ;
      cout.precision(16);
      cout << PreBox << "  GradQ: " << GradientQ(PreBox) <<endl ;
      cout.precision(6);

		}
	}


  cout << "Gradient of Q(Psi|Psi_k) = Q'(" << TestDomain << "|" << Psi_k <<
       ")= " << GradientQ(TestDomain) << endl;

  if (!(0 <= GradientQ(TestDomain))){
		cout << "Gradient of likelihood does not contain zero." << endl;
    cout << "No stationary pt in " << Psi_k << endl << endl ;
	}
  else{

		for (i=1;i<EMiter;i++){

    //--------------------------------------------------------------------
    // First get a copy of Psi_k, then update it for the current iteration
    PrevPsi_k = Psi_k;

			INTERVAL MdPt = IntervalFun(Hull(Mid(Psi_k)));
			cout << "MdPt : " << MdPt << endl;

			REAL LeftTest, RightTest;
			INT MoveLeft, MoveRight;

			INTERVAL Lf = IntervalFun(Hull(Inf(Psi_k)));
			cout << "Lf : " << Lf << endl;
			LeftTest = Inf(Psi_k)-0.5*Diam(Psi_k);
      cout << "LeftTest : " << LeftTest << endl;
			if (LeftTest < Inf(OrigPsi)) LeftTest = Inf(OrigPsi);
      cout << "LeftTest : " << LeftTest << endl;
			INTERVAL TestLf = IntervalFun(Hull(LeftTest));
			cout << "TestLf : " << TestLf << endl;
			if (Inf(TestLf) > Inf(Lf)) MoveLeft=1; else MoveLeft=0;

			INTERVAL Rf = IntervalFun(Hull(Sup(Psi_k)));
			cout << "Rf : " << Rf << endl;
			RightTest = Sup(Psi_k)+0.5*Diam(Psi_k);
      cout << "RightTest : " << RightTest << endl;
      if (RightTest > Sup(OrigPsi)) RightTest=Sup(Psi_k);
      cout << "RightTest : " << RightTest << endl;
			INTERVAL TestRf = IntervalFun(Hull(RightTest));
			cout << "TestRf : " << TestRf << endl;
			if (Inf(TestRf) > Inf(Rf)) MoveRight=1; else MoveRight=0;

      if (MoveLeft==1){
				Psi_k = Hull(LeftTest, Mid(Psi_k));
				cout << "Move Left.  Psi_k: " << Psi_k << endl;
			}
			else if (MoveRight==1){
				Psi_k = Hull(Mid(Psi_k), RightTest);
				cout << "Move Right.  Psi_k: " << Psi_k << endl;
			}
			else{
				if (Inf(Rf)>Inf(Lf)) Psi_k = Hull(Mid(Psi_k),Sup(Psi_k));
				else Psi_k = Hull(Inf(Psi_k),Mid(Psi_k));
				cout << "Bisect.  Psi_k: " << Psi_k << endl;
			}


    // Now intersect the current Psi_k with the original space
    IntersectYN = Intersection(Psi_k,  Psi_k, OrigPsi);
    if (IntersectYN==0) {
      cout << "Note: Psi_k and original Psi_k do not overlap" << endl;
      i = EMiter;
		}
 
    if (Diam(Psi_k) > Diam(PrevPsi_k)){
      Narrow =1;
      // Need to shrink Psi_k.  If it is wider then only one endpoint
      // of the current interval can be inside the previous interval.
      // The other endpoint must be brought back to shrink the interval.
			if ( Inf(Psi_k) <= PrevPsi_k ) 
				Psi_k = Hull( Inf(Psi_k), Inf(Psi_k)+Diam(PrevPsi_k) );
			else
				Psi_k = Hull( Sup(Psi_k)-Diam(PrevPsi_k), Sup(Psi_k) );
		}
    else Narrow = 0;
    if (Narrow==1) 
      cout << "Reduced interval box to previous diameter." << endl;
    cout << "Diam: " << Diam(Psi_k) << endl;
    cout << "Psi_" << i << ": " << Psi_k << endl;

    IVPsi_k = TestDomain;
    IVPsi_k(1) = Psi_k;
    cout << "Q'(Psi_k|Psi_k): " << GradientQ(IVPsi_k) << endl;

    if (Psi_k == PrevPsi_k) {
       i = EMiter;
       cout << endl << "Terminated.  Psi_k = PrevPsi_k." << endl;
	  }

		cout << "==================================" << endl;

	} // End of iteration loop
	} // End of else Gradient
}

