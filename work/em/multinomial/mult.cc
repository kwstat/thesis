// mult.cc
// Multinomial example from DLR

// The vector Psi contains variables: p
// Note, //C is used to denote things that need to be customized/changed

// Need to Code UpdateVars, GradientQ, RealUpper, IntervalUpper, SetTestDomain
// and global constants


#include "UnconstrainedOpt.h"
#include "Functions.h"
//#define DEBUG


// Globals
const INT PROBDIM = 1;  //C
INTERVAL_VECTOR Psi_k(PROBDIM);
INTERVAL C1_k, C2_k;  //C

INTERVAL_VECTOR GradientQ(INTERVAL_VECTOR & Psi)
// Gradient of the likelihood function
{
  static INTERVAL_VECTOR Result(Dimension(Psi));

  INTERVAL P=Psi(1);  //Change variable names and Result vector
  Result(1) = C1_k / P - C2_k / (1-P);

  return (Result);
}



INT GradQNonZero(INTERVAL_VECTOR & Psi)
{
// Check the interval enclosure of the gradient over a box.
// If one direction of the gradient does not contain zero, the box can't 
// contain a local optimum.

  INT NonZero = 0;

  INTERVAL_VECTOR GradQ = GradientQ(Psi);
  INT i, dim=Dimension(Psi);

  for(i=1; i<=dim; i++)
    if ( !(0 <= GradQ(i)) ) NonZero = 1;

  return (NonZero);
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


VOID SetTestDomain(INTERVAL_VECTOR &Psi)
{
  Resize(Psi, PROBDIM);
  cout << "Enter intervals to enclose scalar Psi: ";
  cin >> Psi;
}


INTERVAL IntervalFun(INTERVAL_VECTOR &Psi){

  INTERVAL P=Psi(1);
  REAL q = 0;
  INTERVAL Q = 0;

  Q = C1_k * Log(P) + C2_k * Log(1-P); //C

  Q *= -1.0;

  return (Q);
}

/*
REAL RealUpper (VECTOR &Psi)
// real valued test function
{
  REAL Mu = MuVec(1);
  INT i, lenw = Dimension(W);
  REAL f = 0;
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
    F += Sqr(W(i)-Mu) / (Nu+Sqr(W(i)-Mu_k) );
  F *= (Nu + 1);
  f = Sup(F);

#ifdef DEBUG
  cout << "Real f(Theta):  " << f << endl;
#endif 
  return (f);
}

INTERVAL IntervalUpper (INTERVAL_VECTOR &MuVec)
// interval valued test function
{
  INTERVAL Mu = MuVec(1);
  INT i, lenw = Dimension(W);
  REAL f = 0;
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
    F += Sqr((W(i)-Mu))/(Nu + Sqr((W(i)-Mu_k)));
  F *= (Nu + 1);
  f = Sup(F);

  F = Hull(f);

#ifdef DEBUG
  cout << "Interval f(Theta): " << F << endl;
#endif
  return (F);
}
*/

REAL RealLower(VECTOR &Psi)
// real valued test function
{
  REAL P=Psi(1);
 
  REAL q = 0;
  INTERVAL Q = 0;

  Q = C1_k * Log(P) + C2_k * Log(1-P); //C

  Q *= -1.0;
  q = Inf(Q);

#ifdef DEBUG
  cout << "Real q(Psi):  " << q << endl;
#endif 
  return (q);
}

INTERVAL IntervalLower(INTERVAL_VECTOR &Psi)
// interval valued test function
{

  INTERVAL P=Psi(1);
  REAL q = 0;
  INTERVAL Q = 0;

  Q = C1_k * Log(P) + C2_k * Log(1-P); //C

  Q *= -1.0;
  q = Inf(Q);
  Q = Hull(q);

#ifdef DEBUG
  cout << "Interval q(Psi): " << Q << endl;
#endif
  return (Q);
}

VOID UpdateVars(INTERVAL_VECTOR & Psi_k){
  INTERVAL P_k = Psi_k(1);

  C1_k = 125/(2/P_k +1) + 34;
  C2_k = 18+20;

}

INT main()
{
  INT Iterations, BranchLevels;
  SOLUTIONLIST LowerSolnList, UpperSolnList;
  APPROXIMATIONLIST LowerApproxList, UpperApproxList;
  INTERVAL_VECTOR TestDomain(PROBDIM), PointU, PointL;
  INTERVAL_VECTOR IVPsi_k(PROBDIM);
  REAL LowerBound, UpperBound;

  INT Narrow, i, j, EMiter, IntersectYN;
  INTERVAL_VECTOR PrevPsi_k;

  cout << "Enter max number of EM iterations, EMiter: " ;
  cin >> EMiter;

  cout << "Enter initial interval Psi_k: ";
  cin >> Psi_k;
  UpdateVars(Psi_k);

  INTERVAL_VECTOR OrigPsi = Psi_k;

  SetTestDomain(TestDomain);

  // Added for pre-search, a.k.a. box-grid search
  INT PreSrch;
  REAL BoxWdth, BoxLeft;
  INTERVAL_VECTOR PreBox(1);
  cout << "Enter 0 to skip the presearch" << endl;
  cout << "Enter 1 to do the pre-search." << endl;
  cin >> PreSrch;
  if (PreSrch==1){
    cout << "Enter pre-search box width" << endl;
    cin >> BoxWdth;

    REAL SupPsi = Sup(Psi_k(1));
		BoxLeft=Inf(Psi_k(1));
		do {
      PreBox(1) = Hull(BoxLeft,Sup(Hull(BoxLeft+BoxWdth)));

      BoxLeft = Inf(Hull(BoxLeft+BoxWdth));
		} while (Sup(PreBox(1)) < Sup(Psi_k(1)));

    for (BoxLeft=Inf(Psi_k(1)); BoxLeft<SupPsi; 
           BoxLeft = Inf(Hull(BoxLeft+BoxWdth))) {
      PreBox(1) = Hull(BoxLeft,Sup(Hull(BoxLeft+BoxWdth)));
      if (Sup(PreBox(1)) > SupPsi) 
				PreBox(1) = Hull(Inf(PreBox(1)), SupPsi);
      Psi_k(1)=PreBox(1);
			UpdateVars(Psi_k);
			if (GradQNonZero(PreBox))
        cout << "No   . " ;
			else cout << "Maybe. " ;
      cout.precision(15);
      cout << PreBox << Psi_k << "  GradQ: " << GradientQ(PreBox) <<endl ;
      cout.precision(7);

		}
	}

  cout << "Gradient of Q(Psi|Psi_k) = Q'(" << TestDomain << "|" << Psi_k <<
       ")= " << GradientQ(TestDomain) << endl;


  if (GradQNonZero(TestDomain)){
    cout << "Gradient of likelihood does not contain zero." << endl;
    cout << "No stationary pt in " << Psi_k << endl << endl ;
  }
  else{

    cout << "Computing all global minima of the function" << endl << endl;

    cout << "nit = " << endl; 
    cin >> Iterations;
    cout << "nd  = "; 
    cin >> BranchLevels;
    
    for (i=1;i<EMiter;i++){
/*
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
*/


  //-- Lower --------------------------------------------------------------
  StartUnconstrainedOptimization (LowerSolnList, LowerApproxList,
        Iterations, BranchLevels, 1E-6, 0.2, 0.2,
        LowerBound, UpperBound, TestDomain, RealLower, IntervalLower, 0, 0);

  CleanUpLists (LowerSolnList, LowerApproxList, 1e-6, LowerBound, UpperBound);
    cout << "f min in " << Hull (LowerBound, UpperBound) << endl;
    cout << "LowerSolnList List:" << endl << LowerSolnList << endl;
    if (!IsEmpty(LowerSolnList))
      cout << "Hull of Soln List: " << SolnListHull(LowerSolnList) << endl;
    cout << "LowerApproxList List:" << endl << LowerApproxList << endl;

    if (!IsEmpty(LowerSolnList))
      PointL = SolnListHull(LowerSolnList);
    else
      PointL = Hull(First(LowerApproxList).MinPoint);
    cout << "PointL: " << PointL << endl;


    // Every time through, we need to make sure the lists are empty
    while (!IsEmpty(LowerSolnList)) --LowerSolnList;
    while (!IsEmpty(LowerApproxList)) --LowerApproxList;


    //--------------------------------------------------------------------
    // First get a copy of Psi_k, then update it for the current iteration
    PrevPsi_k = Psi_k;
    

    // Method 1
    // Hull of max point of lower enclosure and endpoint of Psi_k with greater
    // infimum.  Note, intuitively the ">" should be "<", but I'm working
    // with the negative of the Q functions.
    if ( Sup(IntervalFun(Hull(Inf(Psi_k)))) >
            Sup(IntervalFun(Hull(Sup(Psi_k)))) )
      Psi_k = Hull(Sup(Psi_k) ,PointL);
    else
      Psi_k = Hull(Inf(Psi_k), PointL);
    cout << "Left point interval: " << IntervalFun(Hull(Inf(Psi_k))) << endl;
    cout << "Right point interval: " << IntervalFun(Hull(Sup(Psi_k))) << endl;
    cout << "PointL: " << PointL << endl;
    //cout << "Psi_k: " << Psi_k << "   logL: " << LogLike(Psi_k) << endl;


    // Method 2
    // Hull of max point of lower enclosure and max of upper enclosure
    //Mu_k = Hull(PointU(1),PointL(1));


    // Method 3
    // Hull of Midpoint of current interval and max point of lower enclosure
    //Psi_k = Hull(Mid(Psi_k), PointL);



    // Now intersect the current Psi_k with the original space
    IntersectYN = Intersection(Psi_k,  Psi_k, OrigPsi);
    if (IntersectYN==0) {
      cout << "Note: Psi_k and original Psi_k do not overlap" << endl;
      i = EMiter;
    }
 

    //------------------------------------------------------------------
    Narrow =0;
    for (j=1; j<= Dimension(Psi_k); j++){

      if (Diam(Psi_k(j)) > Diam(PrevPsi_k(j))){
        Narrow =1;
        // Need to shrink Psi_k.  If it is wider then only one endpoint
        // of the current interval can be inside the previous interval.
        // The other endpoint must be brought back to shrink the interval.
        if ( Inf(Psi_k(j)) <= PrevPsi_k(j) ) 
          Psi_k(j) = Hull( Inf(Psi_k(j)), Inf(Psi_k(j))+Diam(PrevPsi_k(j)) );
        else
          Psi_k(j) = Hull( Sup(Psi_k(j))-Diam(PrevPsi_k(j)), Sup(Psi_k(j)) );
      }
    }
    UpdateVars(Psi_k);
    if (Narrow==1) cout << "Diameter decreased" << endl;
    else cout << "Diameter same" << endl;
    cout << "Diam: " << Diam(Psi_k) << endl;
    cout.precision(15);
    cout << "Psi_"<<i<<": " << Psi_k << endl;
    cout.precision(6);

//    cout << "Diam PrevPsi_k: " << Diam(PrevPsi_k) << "  Psi_k: " 
//    << Diam(Psi_k) << endl;
//    cout << "   logL: " << LogLike(Mu_k) << endl;

    IVPsi_k = TestDomain;
    cout << "Q'(Psi_k|Psi_k): " << GradientQ(IVPsi_k) << endl;
    cout << endl << "====================================================================" << endl;

  } // End of iteration loop
  } // End of else Gradient
}

