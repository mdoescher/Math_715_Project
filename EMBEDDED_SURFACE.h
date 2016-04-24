//#####################################################################
// Copyright 2014-2016, Michael Doescher, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_SURFACE
//
//#####################################################################
#ifndef __EMBEDDED_SURFACE__
#define __EMBEDDED_SURFACE__

#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include "RANGE_ITERATOR.h"
#include "ELASTICITY_EXAMPLE.h"
#include "REGULAR_HYPERCUBE_MESH.h"
   
namespace PhysBAM{

template<class TV>
class EMBEDDED_SURFACE:public ELASTICITY_EXAMPLE<TV>
{
    enum{d=TV::m};
    typedef typename TV::ELEMENT T;
    typedef ELASTICITY_EXAMPLE<TV> BASE;
    typedef VECTOR<int,d> T_INDEX;
    typedef typename REGULAR_HYPERCUBE_MESH<d>::ELEMENT T_ELEMENT;

protected:
    typedef typename BASE::T_STIFFNESS_MATRIX T_STIFFNESS_MATRIX;

    using BASE::output_directory;
    using BASE::stream_type;
    using BASE::parse_args;
    using BASE::mesh;
    using BASE::X;
    using BASE::X_resting;
    using BASE::stiffness_matrix;
    using BASE::constant_mu;
    using BASE::constant_lambda;
    using BASE::dx;
    using BASE::number_of_newton_iterations;
    using BASE::number_of_CG_iterations;
    using BASE::newton_tolerance;
    using BASE::frame_rate;

  
private:
    T_INDEX size;
    RANGE<TV> domain;
    T_INDEX macroblock_size;
    HASHTABLE<int,int> subdomain_membership; // Element to subdomain
    HASHTABLE<T_INDEX,int> node_hash;        // Auxiliary
    HASHTABLE<int,T_INDEX> node_to_ijk;      // Auxiliary
    HASHTABLE<int> node_is_exterior;
    HASHTABLE<int> node_is_dirichlet;
    HASHTABLE<int> cell_is_full;
//    HASHTABLE<int,T_STIFFNESS_MATRIX> subdomain_stiffness_matrix;

    ARRAY<TV> embedded_vertices;
    ARRAY<VECTOR<int,d> > embedded_faces;
    ARRAY<PAIR<int,TV> > embedding_map;
    ARRAY<VECTOR<T,d> > random_sample_points;
    ARRAY<PAIR<int,TV> > sample_point_embedding_map;
    ARRAY<ARRAY<int> > kinematic_nodes;
    ARRAY<RANGE<TV> > model_ranges;
    int example;
    T rate;
    T perturbation_magnitude;

public:
    EMBEDDED_SURFACE(const STREAM_TYPE stream_type)
        :BASE(stream_type)
    {}

//#####################################################################
public:
    virtual void Write_Output_Files(const int frame) const;
    virtual void Initialize();
    virtual void Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X) const;
    virtual void Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const;
    virtual void Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const;
    virtual void Update_Position_Based_State(const T time);
    void Load_Embedded_Surface();
    void Select_Kinematic_Nodes();
protected:
    virtual void Register_Options();
    virtual void Parse_Options();
//#####################################################################
};

}
#endif

