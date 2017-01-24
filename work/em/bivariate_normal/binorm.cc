// binorm.cc

// Bivariate normal example from EM book and paper by Wu

// Psi = Sigma11, Sigma22, Psi

// Need to Code GradientQ, IntervalFun

#include "UnconstrainedOpt.h"
#include "Functions.h"
#include "gridlist.h"
#include "Utilities.h"


// Globals
const INT PROBDIM = 3;

INTERVAL IntervalFun (INTERVAL_VECTOR &Psi)
{
  INTERVAL S11_k, S22_k, S12_k, Rho_k, Rho2_k, C11, C22, C12;
  S11_k = Psi(1), S22_k=Psi(2), Rho_k=Psi(3);
  S12_k = Rho_k*Sqrt(S11_k*S22_k);
  Rho2_k=Sqr(Rho_k);
  C11 = 20.0+16.0*Sqr(S12_k/S22_k) + 4.0*S11_k*(1-Rho2_k);
  C22 = 20.0+16.0*Sqr(S12_k/S11_k) + 4.0*S22_k*(1-Rho2_k);
  C12 = -16.0*S12_k*(1.0/S22_k + 1.0/S11_k);

  INTERVAL S11, S22, S12, Rho, Det;
  S11=Psi(1), S22=Psi(2), Rho=Psi(3);
  S12=Rho*Sqrt(S11*S22);
  Det = S11*S22*(1-Sqr(Rho));

  INTERVAL Qfun = 0;
  Qfun += -12.0*Log(2.0*Constant::Pi) - 6.0*Log(Det);
  Qfun += -S22/( 2.0*Det ) * C11;
  Qfun += -S11/( 2.0*Det ) * C22;
  Qfun += S12/( Det ) * C12;

  //cout << "Psi: " << Psi << endl;
  //cout << "Q(Psi): " << Qfun << endl;
  return (Qfun);
}

INTERVAL_VECTOR GradientQ(INTERVAL_VECTOR & Psi)
// Gradient of the likelihood function
{

  //cout << endl << endl << "Inside Gradient: " << endl;
  INTERVAL S11_k, S22_k, S12_k, Rho_k, Rho2_k, C11, C22, C12;
  S11_k = Psi(1), S22_k=Psi(2), Rho_k=Psi(3);
  S12_k = Rho_k*Sqrt(S11_k*S22_k);
  Rho2_k=Sqr(Rho_k);
  C11 = 20.0+16.0*Sqr(S12_k/S22_k) + 4.0*S11_k*(1-Rho2_k);
  C22 = 20.0+16.0*Sqr(S12_k/S11_k) + 4.0*S22_k*(1-Rho2_k);
  C12 = 16.0*S12_k*(1.0/S22_k + 1.0/S11_k);

  INTERVAL S11, S22, S12, Rho, Rho2, Omr2, Det;
  S11=Psi(1), S22=Psi(2), Rho=Psi(3);
  S12=Rho*Sqrt(S11*S22);
  Rho2 = Sqr(Rho);
  Omr2 = 1-Rho2;
  Det = S11*S22*(Omr2);
  INTERVAL Opf = Hull(3/2.0);

  static INTERVAL_VECTOR Result(Dimension(Psi));
  INTERVAL D = S11*S22-Sqr(S12);
  Result(1) = -6 + (12+8*Sqr(S12/S22))/(S11*Omr2) -8*Sqr(S12)*(1/S22+1/S11)/(S11*S22*Omr2);
  Result /= S11;
//  Result(1) = -6.0/S11 + C11/(2.0*Sqr(S11)*Omr2)
//            - Rho*C12/(2.0*Power(S11,Opf)*Sqrt(S22)*Omr2);
  Result(2) = -6.0/S22 + C22/(2.0*Sqr(S22)*Omr2)
            - Rho*C12/(2.0*Sqrt(S11)*Power(S22,Opf)*Omr2);
  Result(3) = 12.0*Rho/Omr2 - Rho*C11/(S11*Sqr(Omr2)) - Rho*C22/(S22*Sqr(Omr2))
            + (1+Rho2)*C12/(Sqrt(S11*S22)*Sqr(Omr2));

  //cout << "C11: " << C11 << endl;
  //cout << "-6.0/S11 " << -6.0/S11 << endl;
  //cout << "C11/(2.0*Sqr(S11)*Omr2): " << C11/(2.0*Sqr(S11)*Omr2) << endl;
  //cout << "Rho*C12/(2.0*Power(S11,Opf)*Sqrt(S22)*Omr2): " << Rho*C12/(2.0*Power(S11,Opf)*Sqrt(S22)*Omr2) << endl;
  //cout << "2.0*Power(S11,Opf)*Sqrt(S22)*Omr2: " << 2.0*Power(S11,Opf)*Sqrt(S22)*Omr2 << endl;

  //cout << endl << "C22: " << C22 << endl;
  //cout << "C12: " << C12 << endl;
  //cout << "Result(1): " << Result(1) << endl;
  //cout << "Result(2): " << Result(2) << endl;
  //cout << "Result(3): " << Result(3) << endl;

  //cout << "Psi: " << Psi << endl;
  //cout << "Q'(Psi): " << Result << endl;
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


VOID PrependGrid(GRIDLIST & GridList, INTERVAL_VECTOR &y, INTERVAL &qfun,
  INTERVAL_VECTOR &grad, BOOL &zero){
  GRID_ELEMENT a(y, qfun, grad, zero);
  GridList *= a;
}

VOID QuickSearch(GRIDLIST & GridList){
  First(GridList);
  INTERVAL qfun1, qfun2;
  INTERVAL_VECTOR CurrBox, y1, y2, grad1, grad2;
  INT zero1, zero2, i;
  INT SplitDir, Dim = Dimension(Current(GridList).y);
  REAL MaxDiam, Eps;
  cout << "Max diameter for termination:  " ;
  cin >> Eps;

  INT Done = 0;
  while (!Done){
    First(GridList);
    CurrBox = Current(GridList).y;
    RemoveCurrent(GridList);
    SplitDir = 1;
    for (i=1; i<= Dim; i++)
      if (Diam(CurrBox(i)) > Diam(CurrBox(SplitDir))) SplitDir = i;

      y1 = Lower(CurrBox,SplitDir);
      qfun1 = IntervalFun(y1);
      grad1 = GradientQ(y1);
      zero1 = !(GradQNonZero(y1));

      y2 = Upper(CurrBox,SplitDir);
      qfun2 = IntervalFun(y2);
      grad2 = GradientQ(y2);
      zero2 = !(GradQNonZero(y2));

//      cout << "qfun1: " << qfun1 << "  qfun2 " << qfun2 << endl;
//      cout << "grad1: " << grad1 << "  grad2 " << grad2 << endl;
//      cout << "y1: " << y1 << " y2: " << y2 << endl;

      if (zero1 && !zero2) {
        PrependGrid(GridList,y1,qfun1,grad1,zero1);
      }
      if (zero2 && !zero1){
        PrependGrid(GridList,y2,qfun2,grad2,zero2);
      }
      if (zero1 && zero2){

        if (Inf(qfun2) < Inf(qfun1)) {
          PrependGrid(GridList,y2,qfun2,grad2,zero2);
          PrependGrid(GridList,y1,qfun1,grad1,zero1);
        }
        else {
          PrependGrid(GridList,y1,qfun1,grad1,zero1);
          PrependGrid(GridList,y2,qfun2,grad2,zero2);
        }
      }
       
    
    if (IsEmpty(GridList)) {
      cout << "All boxes eliminated.  Try a bigger initial region." << endl;
      Done=1;
    }
    else {
      First(GridList);
      MaxDiam = 0;
      CurrBox = Current(GridList).y;
      for (i=1; i<= Dim; i++){
        if (Diam(CurrBox(i)) > MaxDiam) 
          MaxDiam = Diam(CurrBox(i));
      }
      if (MaxDiam < Eps) {
        cout << "Eps diameter reached." << endl;
        Done=1;
      }
    }

  } // end of loop
//    cout << "F" << endl;

  return;
}

VOID GridSearch(GRIDLIST & GridList){
  INTERVAL_VECTOR CurrBox, y1, y2, grad1, grad2;
  INTERVAL qenc;
  BOOL zero1, zero2;
  First(GridList);
  INT dim = Dimension(Current(GridList).y);
  INT SplitDir, GridLength, BoxNum;

  // Process the list once
  for (SplitDir = 1; SplitDir <= dim; SplitDir++){
    //cout << "SplitDir: " << SplitDir << endl;
    // Move to start of list and get its length.  Split each box and check.
    First(GridList);
    GridLength = Length(GridList);
    //cout << "GridLength: " << GridLength << endl;
    for (BoxNum = 1; BoxNum <= GridLength; BoxNum++){
      //cout << "BoxNum: " << BoxNum << endl;
      // Deque the box
      CurrBox = Current(GridList).y;
      RemoveCurrent(GridList);
      //cout << "Check first element dequed" << GridList << endl;

      // Split in two pieces, check the enclosure of the gradient
      // Append to the list if necessary
      y1 = Lower(CurrBox,SplitDir);
      grad1 = GradientQ(y1);
      zero1 = !(GradQNonZero(y1));
      if (zero1) {
        qenc = IntervalFun(y1); 
        AppendGrid(GridList, y1, qenc, grad1, zero1);
        //cout << "Appending lower box" << endl;
      }

      y2 = Upper(CurrBox,SplitDir);
      grad2 = GradientQ(y2);
      zero2 = !(GradQNonZero(y2));
      if (zero2) {
        qenc = IntervalFun(y1);
        AppendGrid(GridList, y2, qenc, grad2, zero2);
        //cout << "Appending upper box" << endl;
      }

      //cout << "Check boxes appended: " << GridList << endl;
    }
  }

  return;
}

INT main()
{
  GRIDLIST GridList;  // List of boxes for the grid search
  INTERVAL_VECTOR Psi(PROBDIM);
  INTERVAL_VECTOR Psi_k(PROBDIM);
  INTERVAL Qenc;

  INT i, j;

  cout << "Enter initial interval Psi_k (s11,s22,rho): ";
  cin >> Psi_k;
  Psi = Psi_k;

  IntervalFun(Psi);
  cout << "Gradient of Q(Psi|Psi_k) = " << GradientQ(Psi) << endl;


  if (GradQNonZero(Psi)){
    cout << "Gradient of likelihood does not contain zero." << endl;
    cout << "No stationary pt in " << Psi_k << endl << endl ;
  }
  else{


    //-- Quick Search -------------------------------------------------------
    cout << endl << "Calling Quick Search First" << endl;
    // Prime the list
    Psi_k = Psi;
    INT zero1 = 0;
    // Note, the second Psi and zero1 are just dummy arguments.
    // They are not being used right now.
    AppendGrid(GridList, Psi, Qenc, Psi, zero1);
    QuickSearch(GridList);
    cout.precision(15);
    cout << "After Quick search: " << endl << GridList << endl;
    First(GridList);
    cout << endl << "First element of GridList" << Current(GridList) << endl;


    //-- Bisection Search ----------------------------------------------------

    while (!IsEmpty(GridList)) --GridList;
    // Prime the list
    Psi_k = Psi;
    zero1 = 0;
    // Note, the second Psi and zero1 are just dummy arguments.
    // They are not being used right now.
    AppendGrid(GridList, Psi, Qenc, Psi, zero1);

    cout.precision(16);   
    INT GridSrchIter, TotGridSrchIter=0;
    cout << "Now calling bisection search" << endl;
    cout << "Iterates requested: ";
    cin >> GridSrchIter;
    TotGridSrchIter+=GridSrchIter;
    while (GridSrchIter>0){
      GridSearch(GridList);
      GridSrchIter--;
      if (GridSrchIter==0){
        cout << "GridList: " << GridList << endl;
        cout << "Iterates completed: " << TotGridSrchIter << endl;
        cout << "Add'l iterates: ";
        cin >> GridSrchIter;
        TotGridSrchIter+=GridSrchIter;
      }
    }
    cout << "After GridSearch, here is GridList: " << endl;
    cout << GridList << endl;
    cout.precision(6);


   } // End of else Gradient
}


