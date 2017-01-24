// gridlist.h
// Oct, 1998

#ifndef __GRIDLIST__
#define __GRIDLIST__

#include "IntervalMatrix.h"
#include "Boolean.h"

class GRID_ELEMENT{
  public:
    INTERVAL_VECTOR   y;        // y : the box from the grid search
    INTERVAL          qfun;     // Enclosure of q function
	  INTERVAL_VECTOR   grad;     // grad : enclosure of the gradient over y
    BOOL              zero;     // Does enclosure have zero in some direction?

  GRID_ELEMENT(){};
  GRID_ELEMENT(INTERVAL_VECTOR &,INTERVAL &,INTERVAL_VECTOR &,BOOL &);
  GRID_ELEMENT & operator = (GRID_ELEMENT &);
  friend ostream & operator << (ostream &, GRID_ELEMENT &);
};


#undef LISTOBJECT
#undef LIST
#undef LIST_ELEMENT
#undef LISTCMPFUNC
#define LISTOBJECT   GRIDLISTOBJECT  // What does this do?
#define LIST         GRIDLIST
#define LIST_ELEMENT GRID_ELEMENT
#define LISTCMPFUNC  LCMPFUNC
#include "LinearList.hgen"

// Now ordinary prototypes
VOID AppendGrid(GRIDLIST &, INTERVAL_VECTOR &, INTERVAL &, INTERVAL_VECTOR &,
  BOOL &);


#endif /* __GRIDLIST__ */




