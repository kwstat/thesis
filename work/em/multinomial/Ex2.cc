// Ex2.cc
// Multinomial

#include "Functions.h"   
#include "Utilities.h"

VOID main(){

  INTERVAL PrevT, Theta, Y3 = 0;

  cout << "Initial INTERVAL Theta: " ;
  cin >> Theta;

  INT Converge = 0;
  INT i = 0;
  cout << "i: " << i << " Theta: " << Theta << " Y3: " << Y3 << endl;
  while (Converge == 0){
    i++;
    PrevT = Theta;
//    Y3 = (125 * Theta / 4.0) / ((2 + Theta)/4.0);
    Y3 = 125 / (2/Theta + 1);
    Theta = (34 + Y3)/(72 + Y3);

    cout << "i: " << i << " Theta: " << Theta << " Y3: " << Y3 << endl;
    if (PrevT == Theta) Converge = 1;
	}
}
