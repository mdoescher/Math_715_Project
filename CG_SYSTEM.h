//#####################################################################
// Copyright 2014-2016, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_SYSTEM
//#####################################################################
#ifndef __CG_SYSTEM__
#define __CG_SYSTEM__

//#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>

#include "ELASTICITY_EXAMPLE.h"
#include "CG_VECTOR.h"

namespace PhysBAM{

template<class TV>
class CG_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::ELEMENT>
{
    typedef typename TV::ELEMENT T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;
    
    const ELASTICITY_EXAMPLE<TV>& example;

//    SIMULATION_LAYOUT<T>& layout;
//    const T time;
//    const T dt;

public:
    CG_SYSTEM(const ELASTICITY_EXAMPLE<TV>& example_input)
        :BASE(false,false),example(example_input) {}

    void Multiply(const VECTOR_BASE& x,VECTOR_BASE& y) const
    {
        const ARRAY<TV>& x_array=CG_VECTOR<TV>::Array(x);
        ARRAY<TV>& y_array=CG_VECTOR<TV>::Array(y);
        y_array.Fill(TV());
        example.Add_Force_Differentials(x_array,y_array);
        y_array*=-1.f;
    }

    double Inner_Product(const VECTOR_BASE& x,const VECTOR_BASE& y) const
    {
        const ARRAY<TV>& x_array=CG_VECTOR<TV>::Array(x);
        const ARRAY<TV>& y_array=CG_VECTOR<TV>::Array(y);

        double result=0.;
        for(int i=1;i<=x_array.m;i++)
            result+=TV::Dot_Product(x_array(i),y_array(i));
        return result;
    }

    T Convergence_Norm(const VECTOR_BASE& x) const
    {
        const ARRAY<TV>& x_array=CG_VECTOR<TV>::Array(x);

        T result=0.;
        for(int i=1;i<=x_array.m;i++)
            result=std::max(result,x_array(i).Magnitude());
        return result;
    }

    void Project(VECTOR_BASE& x) const
    {
        ARRAY<TV>& x_array=CG_VECTOR<TV>::Array(x);
        example.Clear_Values_Of_Kinematic_Particles(x_array);
    }

    void Set_Boundary_Conditions(VECTOR_BASE& x) const {}
    
    void Project_Nullspace(VECTOR_BASE& x) const {} // Just a stub (for solids)

//#####################################################################
//#####################################################################
};
}
#endif
