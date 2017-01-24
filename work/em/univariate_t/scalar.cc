// scalar.cc

// EM algorithm, univariate t example of McLachlan, p 95.  
// Scalar case.

#include "UnconstrainedOpt.h"
#include "Functions.h"

// Globals
VECTOR W(4);

INT main()
{
  REAL Eps, Numer, Denom, Mu_k, PrevMu_k;
  INT i, lenw=Dimension(W),Iter,Converge;
  REAL Nu=.05;

  // First define the data
  W(1) = -20;  W(2) = 1; W(3) = 2; W(4) = 3;

  INT Choice;
  cout << " 1. Single starting value.  Print all iterations." << endl;
  cout << " 2. Multiple starts, print final iteration only." << endl;
  cin >> Choice;

  if (Choice==2) {
		cout << "Enter convergence epsilon: ";
		cin >> Eps;

		cout << endl;
		cout << "Convergence epsilon: " << Eps << endl;
		
		Converge = 0;
		Iter = 0;
		REAL StartVal;
		for (StartVal=-20.5;StartVal<=5;StartVal+=.1){
			Mu_k = StartVal;
			Iter = 0;
			Converge = 0;
			
			while(!Converge){
				PrevMu_k = Mu_k;
				Iter++;
				
				// EM step
				Numer = 0;
				Denom = 0;
				for(i=1;i<=lenw;i++)
					Denom += (Nu+1)/(Nu + Sqr(W(i)-Mu_k) );
				for(i=1;i<=lenw;i++)
					Numer += W(i) * (Nu+1)/(Nu + Sqr(W(i)-Mu_k) );
				Mu_k = Numer / Denom;
				
				if (Abs(PrevMu_k - Mu_k) < .0001) Converge = 1;
			}
			cout << "Start: " << StartVal << "     Iterations: " << Iter << "     St. Pt.: " << Mu_k << endl;
		
		} // end StartVal
  }
  else{
		cout << "Enter convergence epsilon: ";
		cin >> Eps;
		cout << "Enter Mu_k: ";
		cin >> Mu_k;
//		cout << "Enter Nu: ";
//		cin >> Nu;
		cout << endl;
		cout << "Convergence epsilon: " << Eps << endl;
    cout << "Mu_0: " << Mu_k << endl;
//    cout << "Nu: " << Nu << endl;
		
		Iter = 0;
		Converge = 0;
		while(!Converge){
			PrevMu_k = Mu_k;
			Iter++;				
			// EM step
			Numer = 0;
			Denom = 0;
			for(i=1;i<=lenw;i++)
				Denom += (Nu+1)/(Nu + Sqr(W(i)-Mu_k) );
			for(i=1;i<=lenw;i++)
				Numer += W(i) * (Nu+1)/(Nu + Sqr(W(i)-Mu_k) );
			Mu_k = Numer / Denom;
			
			if (Abs(PrevMu_k - Mu_k) < .0001) Converge = 1;
			cout << "Iteration: " << Iter << "   Mu_k: " << Mu_k << endl;
		}
			

	}
  return(0);
}
