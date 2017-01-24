// compnorm.cc

#include "IntervalComplex.h"
#include "MiscFunctions.h"    // For the Random() and Rand01() functions.
#include "Functions.h"      // For Cos, Sqr, Constant::Pi
#include "UnconstrainedOpt.h"
#include "Constants.h"        // for Machine::PosInfinity
#include "TestMatrices.h"     // For "matrix" stuff.
#include "Utilities.h"
#include "IntegerVector.h"
#include "Complex.h"
#include "Error.h"

//#define AE

REAL Tolerance;
INTERVAL InvRtPi = Hull(0.5641895835477562,0.5641895835477563);

INTERVAL_COMPLEX Exp(const INTERVAL_COMPLEX& z){
  // z = x+iy
  INTERVAL x = Re(z);
  INTERVAL y = Im(z);
  INTERVAL xRes = Exp(x) * Cos(y);
  INTERVAL yRes = Exp(x) * Sin(y);
  INTERVAL_COMPLEX Result(xRes,yRes);
  return (Result);
}

INTERVAL Gamma2(INT n){
  // This is a special case of the gamma function.
  // Gamma2(n) = Gamma(n+1/2)

  INTERVAL Result;
  INT s,m=2*n-1;
  INTERVAL SqrtPi = Hull(1.77245385090550, 1.77245385090552);
  Result = Hull(1.0);

  while (m>1){
    Result *= m;
    m -= 2;
	}
  Result *= SqrtPi / Power(2, n);
  return(Result);
}

INTERVAL_COMPLEX ErfcAE(const INTERVAL_COMPLEX &z){
  INTERVAL_COMPLEX c_kminus1 = 0.0;
  INTERVAL_COMPLEX c_k = Gamma2(0);  // m=0 term
  INTERVAL_COMPLEX zz = Sqr(z);
  INTERVAL_COMPLEX zzm = 1;
  INTERVAL Converge = Hull(1.0);

  INT m = 0;
  INT Sign = 1;
  while (Inf(Converge) > Tolerance){
    m++;
    Sign *= -1;
    zzm *= zz;
    c_kminus1 = c_k;
    c_k += Sign * Gamma2(m) / zzm;
    Converge = Abs(c_k - c_kminus1); 
    if (m >25) {
      cout << " zzm    " << zzm << endl;
      cout << " Gamma2 " << Gamma2(m) << endl;
      cout << " Sign   " << Sign << endl;
		}

    cout << m << " " << c_k << endl;
	}
  // Now add error part
  REAL Expand = Sup( Gamma2(m+1)/ Abs(zzm * zz));
  INTERVAL RealPart = Hull(Inf(Re(c_k))-Expand, Sup(Re(c_k))+Expand);
  INTERVAL ImagPart = Hull(Inf(Im(c_k))-Expand, Sup(Im(c_k))+Expand);
  INTERVAL_COMPLEX Result(RealPart, ImagPart);

  Result *= Exp(-z*z) / (z * Constant::Pi);
  return (Result);
}


INTERVAL_COMPLEX ContFrac(const INTERVAL_COMPLEX &w){
  // See Jones & Thron, Continued Fractions, p. 316
  // The continued fraction is evaluated with a forward recursive algorithm.
  INTERVAL_COMPLEX p_kminus2 = 1.0;
  INTERVAL_COMPLEX p_kminus1 = 0.0;
  INTERVAL_COMPLEX q_kminus2 = 0.0;
  INTERVAL_COMPLEX q_kminus1 = 1.0;
  INTERVAL_COMPLEX a_k, b_k, c_k, p_k, q_k;
  INTERVAL_COMPLEX c_kminus1 = 0;
  INTERVAL Converge = Hull(1.0);
  //  REAL Tolerance = 1e-10;
  INT MaxIter = 5000;
  // Initializing steps of the continued fraction convergent
  INT k = 0;
  a_k = Hull(1.0);
  b_k = w;
  p_k = b_k * p_kminus1 + a_k * p_kminus2;
  q_k = b_k * q_kminus1 + a_k * q_kminus2;
  c_k = (p_k / q_k);

  while (Inf(Converge) > Tolerance && k < MaxIter){
    k++;
    a_k = Hull(k/2.0);
    b_k = w;
    p_kminus2 = p_kminus1;    
    p_kminus1 = p_k;
    q_kminus2 = q_kminus1;
    q_kminus1 = q_k;
    c_kminus1 = c_k;
    p_k = b_k * p_kminus1 + a_k * p_kminus2;
    q_k = b_k * q_kminus1 + a_k * q_kminus2;
    // Rescale the problem from time to time.  Prevents overflow.
    if (Sup(Abs(p_k)) > 1.0e120){
			p_k /= 1.0e100; p_kminus1 /= 1.0e100;  p_kminus2 /= 1.0e100;
			q_k /= 1.0e100; q_kminus1 /= 1.0e100;  q_kminus2 /= 1.0e100;
    }
    c_k = p_k / q_k;
		Converge = Abs(c_k - c_kminus1);

//  cout << "f_" << k << ": " << c_k << "    " << Converge << endl;
  } 

  // Now expand the box by (c_k - c_k-1)
  REAL Expand = Sup(Converge);  
  INTERVAL RealPart = Hull(Inf(Re(c_k))-Expand, Sup(Re(c_k))+Expand);
  INTERVAL ImagPart = Hull(Inf(Im(c_k))-Expand, Sup(Im(c_k))+Expand);
  INTERVAL_COMPLEX Result(RealPart, ImagPart);

//  cout << "T: " << k+1 << " " ;
  return (Result);
}

INTERVAL_COMPLEX Erfc(const INTERVAL_COMPLEX& w){
  // Complex Error Function

  // Right now this function only handles points which are solely
  // in one quoadrant at a time
  // How can I extend this?

  INTERVAL_COMPLEX z;

  if (Inf(Re(w)) > 0 && Inf(Im(w)) >= 0){  // First quadrant
    z = w;
    return(InvRtPi * Exp(-z*z) * ContFrac(z));
	}
  if (Sup(Re(w)) < 0 && Inf(Im(w)) >=0){ // Second quadrant
    z = -w;
    return(- InvRtPi * Exp(-z*z) * ContFrac(z));
	}
  if (Sup(Re(w)) < 0 && Sup(Im(w))<0){ // Third quadrant
    z = -Conjg(w);
    return(-  Conjg(InvRtPi * Exp(-z*z) * ContFrac(z))  );
	}
  if (Inf(Re(w)) >0 && Inf(Im(w)) <=0){  // Fourth quadrant
    z = Conjg(w);
    return(  Conjg(InvRtPi * Exp(-z*z) * ContFrac(z))  );
	}
}

INTERVAL_COMPLEX w(const INTERVAL_COMPLEX &z){
  // w(z)  Abramowitz&Stegun p 325.

  INTERVAL_COMPLEX Result = Exp(-z*z) * Erfc(-I*z);
  return(Result);
}

INTERVAL_COMPLEX P(const INTERVAL_COMPLEX &z){
  INTERVAL_COMPLEX Result = 0.5 * Exp(-z*z) * w(-I*z/Sqrt(2.0));
  return(Result);
}

INT main(){
  INT Choice;
  cout.precision(17);

  cout << "1.  Evaluate for a single point" << endl;
  cout << "2.  Reproduce Abram./Stegun table" << endl;
  cin >> Choice;

  cout << "Enter tolerance: " ;
  cin >> Tolerance;
  cout << endl;

  if (Choice==1){
		INTERVAL xr,xi,yr,yi;
		INTERVAL_COMPLEX x,y,z;
		REAL  rl, ru, il, iu;
		cout << "Input real_l, real_u, imag_l, imag_u: ";
		cin >> rl >> ru >> il >> iu;
		cout << endl;
		INTERVAL zRe(rl,ru), zIm(il,iu);
		z = INTERVAL_COMPLEX(zRe,zIm);
		cout << zRe << " " << zIm << " " << z << endl;
		cout << "w(z): " << w(z) << endl;
		cout << "P(z): " << P(z) << endl;
	}

  if (Choice==2){

		INT Space = 0;
    INT xind, yind;
		REAL x,y;
		for (xind=1;xind<=39; xind++){
      x = xind/10.0;
			cout << endl;
			cout << "x: " << x << endl;

      for (yind=1;yind<=30; yind++){
				y = yind/10.0;

				INTERVAL_COMPLEX z(x,y);
				//		  cout << "w("<<z<<") " << w(z) << endl;
				cout << y << "  " << w(z) << endl;
				Space++;
				if (Space == 5){
					Space = 0;
					cout << endl;
				}
			}
		}
	}

#ifdef AE
  INT m;
  cout << "A " << endl;
  for (m=0; m<= 10; m++){
    cout << "Gamma2(" << m << ") " << Gamma2(m) << endl;
	}
  cout << "----------------" << endl;
  
  INTERVAL_COMPLEX z;
  REAL  rl, ru, il, iu;
  cout << "Input rl, ru, il, iu: ";
  cin >> rl >> ru >> il >> iu;
  cout << endl;
  INTERVAL zRe(rl,ru), zIm(il,iu);
  z = INTERVAL_COMPLEX(zRe,zIm);
  cout << "ErfcAE("<<z<<") " << ErfcAE(z) << endl;
#endif // AE
}





