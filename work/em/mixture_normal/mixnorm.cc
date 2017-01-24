// mixnorm.cc

// Mixture of two normal dist'ns from EM book

// Psi = Mu1, Mu2, Sig, P

// Need to Code GradientQ, IntervalFun

#include "UnconstrainedOpt.h"
#include "Functions.h"
#include "gridlist.h"
#include "Utilities.h"


// Globals
const INT PROBDIM = 4;
const INT DATADIM = 20;
INTERVAL_VECTOR W(DATADIM);

INTERVAL IntervalFun (INTERVAL_VECTOR &Psi)
{

  INTERVAL Qfun = 0;

  //cout << "Psi: " << Psi << endl;
  //cout << "Q(Psi): " << Qfun << endl;
  return (Qfun);
}

INTERVAL_VECTOR GradientQ(INTERVAL_VECTOR & Psi)
// Gradient of the likelihood function
{
  static INTERVAL_VECTOR Result(PROBDIM);

  //cout << endl << endl << "Inside Gradient: " << endl;
  static INTERVAL_VECTOR Z1(DATADIM), Z2(DATADIM);
  INT i,j;
  INTERVAL Mu1 = Psi(1), Mu2=Psi(2), Sig=Psi(3), P=Psi(4);
  INTERVAL E, Sig2 = Sqr(Sig), Sig3=Sig*Sig2;

  for(j=1;j<=DATADIM;j++){
    E = Exp((-Sqr(W(j)-Mu2)+Sqr(W(j)-Mu1))/(2*Sig2) );
    Z1(j) = P/(P+(1-P)*E);
    E = Exp((-Sqr(W(j)-Mu1)+Sqr(W(j)-Mu2))/(2*Sig2));
    Z2(j) = (1-P)/(P*E+1-P);
	}

  Result(1)=0;Result(2)=0;Result(3)=0;Result(4)=0;
  for(j=1;j<=DATADIM;j++){
    Result(1) += Z1(j)*(W(j)-Mu1);
    Result(2) += Z2(j)*(W(j)-Mu2);
    Result(3) += Z1(j)*(-1/Sig +Sqr(W(j)-Mu1)/Sig3)
                +Z2(j)*(-1/Sig +Sqr(W(j)-Mu2)/Sig3);
    Result(4) += Z1(j)/P -Z2(j)/(1-P);
	}
  Result(1)/=Sig2;
  Result(2)/=Sig2;

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

// Hardcode data for now.
W(1)=1.30;W(2)=-1.16;
W(3)= 0.31;W(4)= 0.88;W(5)= 1.40;W(6)= 2.86;W(7)= 1.30;W(8)= 2.77;
W(9)= 1.79;W(10)= 2.22;W(11)= 2.41;W(12)= 0.94;W(13)= 1.75;W(14)= 2.68;
W(15)= 2.28;W(16)= 2.51;W(17)= 2.24;W(18)= 2.75;W(19)= 1.82;W(20)= 2.10;

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


