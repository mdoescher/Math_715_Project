//#####################################################################
// Copyright 2014-2016, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTICITY_DRIVER
//#####################################################################
#ifndef __ELASTICITY_DRIVER__
#define __ELASTICITY_DRIVER__

#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include "ELASTICITY_EXAMPLE.h"

namespace PhysBAM{

template<class TV>
class ELASTICITY_DRIVER:public DRIVER<TV>
{
    enum{d=TV::m};
    typedef typename TV::ELEMENT T;
    typedef DRIVER<TV> BASE;
    ELASTICITY_EXAMPLE<TV>& example;

public:
    ELASTICITY_DRIVER(ELASTICITY_EXAMPLE<TV>& example_input)
        :BASE(example_input),example(example_input)
    {}

//#####################################################################
    virtual void Initialize();
    virtual void Advance_To_Target_Time(const T target_time);
//#####################################################################
};
}
#endif
