// gridlist.cc
// Oct, 1998

#include "gridlist.h"
#include "Constants.h"
#include "Utilities.h"
#include "LinearList.Cgen"

#undef LISTOBJECT
#undef LIST
#undef LIST_ELEMENT
#undef LISTCMPFUNC

// Constructor
GRID_ELEMENT::GRID_ELEMENT(INTERVAL_VECTOR &ybox, INTERVAL &qenc,
  INTERVAL_VECTOR & gbox, BOOL & zeroval)
{
  y = ybox;
  qfun = qenc;
  grad = gbox;
  zero = zeroval;
}

// Assignment operator
GRID_ELEMENT & GRID_ELEMENT::operator = (GRID_ELEMENT &v)
{
  y = v.y;
  qfun = v.qfun;
  grad = v.grad;
  zero = v.zero;
}

// Out-stream operator
ostream & operator << (ostream & o, GRID_ELEMENT &v)
{
  o << endl;
  o << "y: " << v.y << endl;
  o << "qfun: " << v.qfun << endl;
  o << "grad: " << v.grad << endl;
  return o;
}

// Append element to end of list
VOID AppendGrid(GRIDLIST &list, INTERVAL_VECTOR &y, INTERVAL &qfun, 
  INTERVAL_VECTOR &grad, BOOL &zero)
{
  // First define an element, then add it to the list
  GRID_ELEMENT a(y, qfun, grad, zero);
  list += a;
}


