// norm.cc
// Bivariate normal example from EM book and paper by Wu
// Psi = Sigma11, Sigma22, Psi

#include "UnconstrainedOpt.h"
#include "Functions.h"
//#define DEBUG


// Globals
INTERVAL_VECTOR Psi_k(3);
INTERVAL C11, C22, C12;

/*
// Evaluate the log likelihood for Mu
INTERVAL LogLike(INTERVAL & Mu)
{
  INTERVAL Result = 0;
  return (Result);
}
*/

/*
INTERVAL_VECTOR GradientLike(INTERVAL_VECTOR & Psi)
{

}
*/

INTERVAL_VECTOR GradientQ(INTERVAL_VECTOR & Psi)
// Gradient of the likelihood function
{
  INTERVAL S11=Psi(1), S22=Psi(2), Rho=Psi(3);
  INTERVAL S12=Rho*Sqrt(S11*S22);
  INTERVAL Rho2=Sqr(Rho);
  INTERVAL Omr2 = 1-Rho2;
  INTERVAL Det = S11*S22*(1-Rho2);
//  cout << "Det " << Det << endl;


  static INTERVAL_VECTOR Result(Dimension(Psi));
//  cout << "Result: " << Result << endl;

  Result(1) = -6.0/S11 + C11/(2.0*Sqr(S11)*Omr2)
            + Rho*C12/(2.0*Power(S11,3/2)*Sqrt(S22)*Omr2);
  Result(2) = -6.0/S22 + C22/(2.0*Sqr(S22)*Omr2)
            + Rho*C12/(2.0*Sqrt(S11)*Power(S22,3/2)*Omr2);
  Result(3) = 12.0*Rho/Omr2 - Rho*C11/(S11*Sqr(Omr2)) - Rho*C22/(S22*Sqr(Omr2))
            - (1+Rho2)*C12/(Sqrt(S11*S22)*Sqr(Omr2));

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


VOID SetTestDomain (INTERVAL_VECTOR &Psi)
{
  Resize(Psi, 3);
  cout << "Enter intervals to enclose scalar Psi: ";
	cin >> Psi;
}

/*
INTERVAL IntervalFun (INTERVAL &Psi){
  INT i, lenw = Dimension(W);
  INTERVAL F = 0;
  for (i=1;i<=lenw;i++)
    F += Sqr(W(i)-Mu) / (Nu+Sqr(W(i)-Mu_k) );
  F *= (Nu + 1);
  return(F);
}
*/
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

REAL RealLower (VECTOR &Psi)
// real valued test function
{
  REAL S11=Psi(1), S22=Psi(2), Rho=Psi(3);
  INTERVAL S12=Rho*Sqrt(S11*S22);
  INTERVAL Rho2 = Sqr(Rho);
  INTERVAL Det = S11*S22*(1-Rho2);
 
  REAL q = 0;
  INTERVAL Q = 0;
  if (Inf(Det)==0) Det = Hull(.01, Sup(Det));  // HACK HACK HACK
#ifdef DEBUG
  cout << "RealLower" << endl;
  cout << "Psi: " << Psi << endl;
  cout << "Det: " <<  Det << endl;
#endif
  Q += -12.0*Log(2.0*Constant::Pi) - 6.0*Log(Det);
  Q += -S22/( 2.0*Det ) * C22;
  Q += -S11/( 2.0*Det ) * C11;
  Q += -S12/( Det ) * C12;
  Q *= -1.0;  // Software is set up for minimization, not maximization

// Need to calculate C11, etc.
  q = Inf(Q);

#ifdef DEBUG
  cout << "Real q(Psi):  " << q << endl;
#endif 
  return (q);
}

INTERVAL IntervalLower (INTERVAL_VECTOR &Psi)
// interval valued test function
{
  INTERVAL S11=Psi(1), S22=Psi(2), Rho=Psi(3);
  INTERVAL S12=Rho*Sqrt(S11*S22);
  INTERVAL Rho2=Sqr(Rho);
  INTERVAL Det = S11*S22*(1-Rho2);

  REAL q = 0;
  INTERVAL Q = 0;
  if (Inf(Det)==0) Det = Hull(.01, Sup(Det));  // HACK HACK HACK
#ifdef DEBUG
  cout << "IntervalLower" << endl;
  cout << "Psi: " << Psi << endl;
  cout << "Det: " <<  Det << endl;
#endif
  Q += -12.0*Log(2.0*Constant::Pi) - 6.0*Log(Det);
  Q += -S22/( 2.0*Det ) * C22;
  Q += -S11/( 2.0*Det ) * C11;
  Q += -S12/( Det ) * C12;
  Q *= -1.0;

  q = Inf(Q);
  Q = Hull(q);

#ifdef DEBUG
  cout << "Interval q(Psi): " << Q << endl;
#endif
  return (Q);
}

VOID UpdateVars(INTERVAL_VECTOR & Psi_k){

  INTERVAL S11_k=Psi_k(1), S22_k=Psi_k(2), S12_k=Psi_k(3);
  C11 = 20.0+16.0*Sqr(S12_k/S22_k) + 4.0*S22_k*(1-Sqr(S12_k)/(S11_k * S22_k));
  C22 = 20.0+16.0*Sqr(S12_k/S11_k) + 4.0*S11_k*(1-Sqr(S12_k)/(S11_k * S22_k));
  C12 = 16.0*S12_k*(1.0/S22_k + 1.0/S11_k);

}

INT main()
{
  INT Iterations, BranchLevels;
  SOLUTIONLIST LowerSolnList, UpperSolnList;
  APPROXIMATIONLIST LowerApproxList, UpperApproxList;
  INTERVAL_VECTOR TestDomain(3), PointU, PointL;
  INTERVAL_VECTOR IVPsi_k(3);
  REAL LowerBound, UpperBound;

  INT Narrow, i, j, EMiter, IntersectYN;
  INTERVAL_VECTOR PrevPsi_k;

  cout << "Enter number of EM iterations, EMiter: " ;
  cin >> EMiter;

  cout << "Enter initial interval Psi_k: ";
  cin >> Psi_k;
  UpdateVars(Psi_k);

  INTERVAL_VECTOR OrigPsi = Psi_k;

	SetTestDomain(TestDomain);

//  cout << endl << "Gradient of LogL(Psi_0) = LogL'(" << Psi_k << "): " <<
//       GradientLike(Psi_k) << endl;

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
//		cout << "UpperSolnList List:" << endl << UpperSolnList << endl;
		if (!IsEmpty(UpperSolnList))
			cout << "Hull of Soln List: " << SolnListHull(UpperSolnList) << endl;
//		cout << "UpperApproxList List:" << endl << UpperApproxList << endl;

		if (!IsEmpty(UpperSolnList))
			PointU = SolnListHull(UpperSolnList);
		else
			PointU = Hull(First(UpperApproxList).MinPoint);
		cout << "PointU: " << PointU << endl;

    // Every time through, we need to make sure the lists are empty
    while (!IsEmpty(UpperSolnList)) --UpperSolnList;
    while (!IsEmpty(UpperApproxList)) --UpperApproxList;
*/


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

		// Every time through, we need to make sure the lists are empty
		while (!IsEmpty(LowerSolnList)) --LowerSolnList;
		while (!IsEmpty(LowerApproxList)) --LowerApproxList;



    //--------------------------------------------------------------------
    // First get a copy of Psi_k, then update it for the current iteration
    PrevPsi_k = Psi_k;
    
/*
    // Hull of max point of lower enclosure and endpoint of Psi_k with greater
    // infimum.  Note, intuitively the ">" should be "<", but I'm working
    // with the negative of the Q functions.
    if ( Sup(IntervalFun(Hull(Inf(Mu_k)))) > Sup(IntervalFun(Hull(Sup(Mu_k)))) )
			Mu_k = Hull(Sup(Mu_k) ,PointL(1));
		else
			Mu_k = Hull(Inf(Mu_k), PointL(1));
    cout << "Left point interval: " << IntervalFun(Hull(Inf(Mu_k))) << endl;
    cout << "Right point interval: " << IntervalFun(Hull(Sup(Mu_k))) << endl;
    cout << "PointL(1): " << PointL(1) << endl;
		cout << "Mu_k: " << Mu_k << "   logL: " << LogLike(Mu_k) << endl;
*/

    // Hull of max point of lower enclosure and max of upper enclosure
    //Mu_k = Hull(PointU(1),PointL(1));

    // Hull of Midpoint of current interval and max point of lower enclosure
    Psi_k = Hull(Mid(Psi_k), PointL);


    //-------------------------------------------------------------
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
    cout << Narrow << " Reduced Diam: " << Diam(Psi_k) << endl;
    cout << " Psi_k: " << Psi_k << endl;

//    cout << "Diam PrevPsi_k: " << Diam(PrevPsi_k) << "  Psi_k: " 
//    << Diam(Psi_k) << endl;
//    cout << "   logL: " << LogLike(Mu_k) << endl;

    IVPsi_k = TestDomain;
    cout << "Q'(Psi_k|Psi_k): " << GradientQ(IVPsi_k) << endl;
		cout << endl << "====================================================================" << endl;

	} // End of iteration loop
	} // End of else Gradient
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
/*
Gradient of Q parameterized by S11, S22, S12

  Result(1) = -6.0* S22/Det    +Sqr(S22/Det)*C11/2.0
              +Sqr(S12/Det)*C22/2.0   + S12*S22*C12/Sqr(Det);
  Result(2) = -6.0*S11/Det   +Sqr(S12/Det)*C11/2.0
              +Sqr(S11/Det) + S12*S11/Sqr(Det)*C11;
  Result(3) = 12.0*S12/Det  -S22*S12*C11/Sqr(Det)
              -S11*S12*C22/Sqr(Det)  -(S11*S22+Sqr(S12))*C12/Sqr(Det);



*/
