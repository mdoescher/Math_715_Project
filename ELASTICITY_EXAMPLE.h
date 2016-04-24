//#####################################################################
// Copyright 2014-2016, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTICITY_EXAMPLE
//#####################################################################
#ifndef __ELASTICITY_EXAMPLE__
#define __ELASTICITY_EXAMPLE__

#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>

namespace PhysBAM{

template<int d> class REGULAR_HYPERCUBE_MESH;

template<class TV>
class ELASTICITY_EXAMPLE:public EXAMPLE<TV>
{
    enum{d=TV::m};
    typedef typename TV::ELEMENT T;
    typedef EXAMPLE<TV> BASE;

protected:
    typedef HASHTABLE<VECTOR<int,2>,MATRIX<T,d> > T_STIFFNESS_MATRIX;

    using BASE::output_directory;
    using BASE::stream_type;
    using BASE::parse_args;
    using BASE::frame_rate;

public:
    REGULAR_HYPERCUBE_MESH<d> mesh;
    ARRAY<TV> X;
    ARRAY<TV> X_resting;

    T constant_mu;
    T constant_lambda;
    T dx;
    int number_of_newton_iterations;
    int number_of_CG_iterations;
    T newton_tolerance;

    T_STIFFNESS_MATRIX stiffness_matrix;

    ELASTICITY_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type)
    {}

//#####################################################################
public:
    virtual void Initialize();
    virtual void Write_Output_Files(const int frame) const;
    virtual void Update_Position_Based_State(const T time);
    virtual void Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const;
    virtual void Add_Force_Differentials(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dforce) const;
    virtual void Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X) const = 0;
    virtual void Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const = 0;
    virtual void Update_Collisions(const T time){}
protected:
    virtual void Register_Options();
    virtual void Parse_Options();
//#####################################################################
};
}
#endif
