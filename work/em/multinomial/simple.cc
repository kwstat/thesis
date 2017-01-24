// simple.cc
// The multinomial example of Dempster, Laird, and Rubin.
// (1)  The scalar algorithm of DLR.
// (2)  The same algorithm, using intervals instead of scalars.
// (2a) The original formulation
// (2b) A formulation with less dependency.

#include <iostream.h>
#include <MiscFunctions.h>  // For the Random() and Rand01() functions.
#include <Functions.h>      // For the Exp function.
#include <math.h>           // For the pow,cos,  functions.
#include <stdlib.h>
#include <Vector.h>         // Definitions of VECTOR type
#include <Interval.h>

void Multinomial(void)
{
  INTERVAL p, x2;
  int Converge;
  double Eps = 0.0000001;
  double Pi, dx2;
  double PrevPi;
  int Iteration;

  // First the non-interval case ----------------------------------------

  Pi = 0;
  cout << "Multinomial example from DLR" << endl;
  cout << "First the scalar case.  Input parameter: p" << endl;
  while (1) {
    cout << "Enter a number greater than 1 to exit scalar case." << endl;
    cout << "Enter pi(0):  ";
    cin >> Pi;
    if (Pi < 0 || Pi > 1) break;
    cout << "Epsilon:    " << Eps << endl;
    cout << "Initial pi: " << Pi << endl; 
    cout << "Iter     pi           x2       " << endl;

    Iteration = 0;
    Converge = 0;
    while (Converge != 1){
      ++Iteration;
      PrevPi = Pi;
      dx2 = 125.0 * Pi / (2.0 + Pi);
      Pi = 1.0 - 38.0 / (dx2 + 72.0);

      cout << Iteration << "   " << Pi << "   " << dx2 << endl;
      // Now check for convergence    
      if ( fabs(PrevPi-Pi) < Eps ) Converge = 1;
    } 
    cout << "Converges when |PrevPi-Pi| < Eps." << endl;
  }

  //-------------------------------------------------------------------------
  // Now the interval case.  The first formulation (with dependency)
  cout << endl << endl;
  cout << "First interval formulation (with dependency)" << endl;
  p = Hull(0.0,1.0);
  cout << "Epsilon:   " << Eps << endl;
  cout << "Initial p: " << p << endl;
 
  cout << "Iter     p           x2       " << endl;

  Iteration = 0;
  Converge = 0;
  while (Converge != 1){
    ++Iteration;
    x2 = 125.0 * p / (2.0 + p);
    p = 1.0 - 38.0 / (x2 + 72.0);

    cout << Iteration << " " << p << " " << x2 << endl;

    // Now check for convergence
    
    if ( Diam(p) < Eps ) Converge = 1;
  }
  cout << "Converges when Diam(p) < Eps." << endl;


  //-------------------------------------------------------------------------
  // Interval case, different formulation with p = (0,1], no dependency
  cout << endl << endl;
  cout << "Second interval formulation (with dependency)" << endl;
  p = Hull(Succ(0.0), 1.0);
  cout << "Epsilon:   " << Eps << endl;
  cout << "Initial p: " << p << endl;
 
  cout << "Iter     p           x2       " << endl;

  Iteration = 0;
  Converge = 0;
  while (Converge != 1){
    ++Iteration;
    x2 = 125.0 / (2 / p + 1);
    p = 1.0 - 38.0 / (x2 + 72.0);

    cout << Iteration << " " << p << " " << x2 << endl;

    // Now check for convergence
    
    if ( Diam(p) < Eps ) Converge = 1;
  } 
  cout << "Converges when Diam(p) < Eps." << endl;
  return;
}

int main(void){
  Multinomial();
  return(0);
}
