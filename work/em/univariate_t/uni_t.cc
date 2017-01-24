// uni_t.cc
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


// Evaluate the log likelihood for Psi
INTERVAL LogLike(INTERVAL & Psi)
{
  INTERVAL Mu = Psi;
  INTERVAL Result = 0;
  INT i, lenw=Dimension(W);
  
  for(i=1;i<=lenw;i++){
    Result += Log(1+20*Sqr(W(i)-Mu));
	}
  Result *= -1.0;
  return (Result);
}

INTERVAL GradientLike(INTERVAL & Psi)
// Gradient of the likelihood function
{
  INTERVAL Mu = Psi;
  INTERVAL Result = 0;
  INT i, lenw=Dimension(W);

  for(i=1;i<=lenw;i++){
		Result += 2*(W(i)-Mu)/(Nu+ Sqr(W(i)-Mu) ) ;
	}
  return (Result);
}

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

INTERVAL_VECTOR SolnListHull(SOLUTIONLIST & S)
{
	First(S);
	INTERVAL_VECTOR Box = Current(S).Box;
	while (!Finished(S)){
		Box = Hull(Box,Current(S).Box);
		Next(S);
	}
	
  return (Box);
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
  F *= (Nu + 1);
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
  F *= (Nu + 1);
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
  F *= (Nu + 1);
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

  INT Iterations, BranchLevels;
  SOLUTIONLIST LowerSolnList, UpperSolnList;
  APPROXIMATIONLIST LowerApproxList, UpperApproxList;
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


  cout << endl << "Gradient of LogL(Psi_0) = LogL'(" << Psi_k << "): " <<
       GradientLike(Psi_k) << endl;
  cout << "Gradient of Q(Psi|Psi_k) = Q'(" << TestDomain << "|" << Psi_k <<
       ")= " << GradientQ(TestDomain) << endl;

  if (!(0 <= GradientQ(TestDomain))){
		cout << "Gradient of likelihood does not contain zero." << endl;
    cout << "No stationary pt in " << Psi_k << endl << endl ;
	}
  else{

    cout << "Enter method." << endl;
    cout << "1 = Max of lower, endpoint" << endl;
    cout << "2 = Max of lower, max of upper (not working)" << endl;
                cout << "3 = Mid point of current, max of lower" << endl;
    cin >> Method;

		cout << "Computing all global minima of the function" << endl << endl;

		cout << "nit = " << endl; 
		cin >> Iterations;
		cout << "nd  = "; 
		cin >> BranchLevels;
		
		for (i=1;i<EMiter;i++){

/* if (Method==2){
			//---- Upper ---------------------------------------------------
			cout << "Upper -----------------------" << endl;
			StartUnconstrainedOptimization (UpperSolnList, UpperApproxList,
					Iterations, BranchLevels, 1E-6, 0.2, 0.2,
				  LowerBound, UpperBound, TestDomain,
				  RealUpper, IntervalUpper, 0, 0);

		CleanUpLists (UpperSolnList, UpperApproxList, 1e-6, LowerBound,UpperBound);
    cout << "f min in " << Hull (LowerBound, UpperBound) << endl;
    //cout << "UpperSolnList List:" << endl << UpperSolnList << endl;
		if (!IsEmpty(UpperSolnList))
			cout << "Hull of Soln List: " << SolnListHull(UpperSolnList) << endl;
    //cout << "UpperApproxList List:" << endl << UpperApproxList << endl;

		if (!IsEmpty(UpperSolnList))
			PointU = SolnListHull(UpperSolnList);
		else
			PointU = Hull(First(UpperApproxList).MinPoint);
		cout << "PointU: " << PointU << endl;

    // Every time through, we need to make sure the lists are empty
    while (!IsEmpty(UpperSolnList)) --UpperSolnList;
    while (!IsEmpty(UpperApproxList)) --UpperApproxList;
 } */


		//--------------------------------------------------------------
		cout << "Lower -----------------------" << endl;
		StartUnconstrainedOptimization (LowerSolnList, LowerApproxList,
				  Iterations, BranchLevels, 1E-6, 0.2, 0.2,
				  LowerBound, UpperBound, TestDomain,
				  RealLower, IntervalLower, 0, 0);

  CleanUpLists (LowerSolnList, LowerApproxList, 1e-6, LowerBound, UpperBound);
    cout << "f min in " << Hull (LowerBound, UpperBound) << endl;
//		cout << "LowerSolnList List:" << endl << LowerSolnList << endl;
		if (!IsEmpty(LowerSolnList))
			cout << "Hull of Soln List: " << SolnListHull(LowerSolnList) << endl;
//		cout << "LowerApproxList List:" << endl << LowerApproxList << endl;

		if (!IsEmpty(LowerSolnList))
			PointL = SolnListHull(LowerSolnList);
		else
			PointL = Hull(First(LowerApproxList).MinPoint);
		cout << "PointL: " << PointL << endl;



    //--------------------------------------------------------------------
    // First get a copy of Psi_k, then update it for the current iteration
    PrevPsi_k = Psi_k;
    
		// Every time through, we need to make sure the lists are empty
		while (!IsEmpty(LowerSolnList)) --LowerSolnList;
		while (!IsEmpty(LowerApproxList)) --LowerApproxList;

    if (Method==1) {
			// Hull of max point of lower enclosure and endpoint of Psi_k with greater
			// infimum.  Note, intuitively the ">" should be "<", but I'm working
			// with the negative of the Q functions.
			if ( Sup(IntervalFun(Hull(Inf(Psi_k)))) > 
                        Sup(IntervalFun(Hull(Sup(Psi_k)))) )
				Psi_k = Hull(Sup(Psi_k) ,PointL(1));
			else
				Psi_k = Hull(Inf(Psi_k), PointL(1));
			cout << "Left point interval: " << IntervalFun(Hull(Inf(Psi_k))) << endl;
			cout << "Right point interval: " << IntervalFun(Hull(Sup(Psi_k))) << endl;
			cout << "PointL(1): " << PointL(1) << endl;
	 		cout << "Psi_k: " << Psi_k << "   logL: " << LogLike(Psi_k) << endl;
		}

    /* if (Method==2){
      // Hull of max point of lower enclosure and max of upper enclosure
      //Psi_k = Hull(PointU(1),PointL(1));
    */

    if (Method==3){
      // Hull of Midpoint of current interval and max point of lower enclosure
      Psi_k = Hull(Mid(Psi_k), PointL(1));
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

    cout << "   logL: " << LogLike(Psi_k) << endl;
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

