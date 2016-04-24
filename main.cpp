// #include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
// #include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
// #include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
// #include "RANGE_ITERATOR.h"

// #include "REGULAR_HYPERCUBE_MESH.h"

#include "EMBEDDED_SURFACE.h"
#include "ELASTICITY_DRIVER.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef float T;
    typedef float RW;
    STREAM_TYPE stream_type((RW()));
    enum {d=3};
    typedef VECTOR<T,d> TV;
    LOG::Initialize_Logging();
    EMBEDDED_SURFACE<TV> example(stream_type);
    example.Parse(argc,argv);
    ELASTICITY_DRIVER<TV> driver(example);
    driver.Execute_Main_Program();

    LOG::Finish_Logging();
    return 0;
}
