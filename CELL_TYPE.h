#ifndef __CELL_TYPE__
#define __CELL_TYPE__

namespace PhysBAM
{

typedef enum {INTERIOR_CELL_TYPE=1,
              DIRICHLET_CELL_TYPE=2,
              EXTERIOR_CELL_TYPE=3,
              BOUNDARY_CELL_TYPE=4,
              PATHOLOGICAL_CELL_TYPE=5} CELL_TYPE;


}

#endif
