//#####################################################################
// Copyright 2014-2016, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTICITY_DRIVER
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include "ELASTICITY_DRIVER.h"
#include "REGULAR_HYPERCUBE_MESH.h"
#include "CG_VECTOR.h"
#include "CG_SYSTEM.h"
using namespace PhysBAM;
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void ELASTICITY_DRIVER<TV>::
Initialize()
{
    example.Parse_Late_Options();
    BASE::Initialize();
    example.Initialize();
    example.Write_Output_Files(0);
}
//#####################################################################
// Function Advance_To_Target_Time
//#####################################################################
template<class TV> void ELASTICITY_DRIVER<TV>::
Advance_To_Target_Time(const T target_time){
    LOG::SCOPE scope("{Advance_To_Target_Time");
   
        example.Update_Collisions(target_time);
        example.Set_Kinematic_Positions(target_time,example.X);

    for (int newton_iteration=1;newton_iteration<=example.number_of_newton_iterations;newton_iteration++) {
        ARRAY<TV> dX(example.X.Size()),forces(example.X.Size()),q(example.X.Size()),s(example.X.Size()),r(example.X.Size());

        // Construct right hand side
        example.Update_Position_Based_State(target_time);
        example.Add_Elastic_Forces(example.X,forces);
        example.Clear_Values_Of_Kinematic_Particles(forces);
        
        // Krylov solver wrappers
        CG_VECTOR<TV> cg_x(dX),cg_b(forces),cg_q(q),cg_s(s),cg_r(r);
        CG_SYSTEM<TV> cg_system(example);
        
        T norm=cg_system.Convergence_Norm(cg_b);
        LOG::cout<<newton_iteration<<" Forces norm "<<norm<<std::endl;
        if (norm <=example.newton_tolerance) break;

        // Generate Conjugate Gradients solver object
        CONJUGATE_GRADIENT<T> cg;
        cg.print_residuals=false;
        cg.print_diagnostics=true;

        // Solve linear system using CG
        cg.Solve(cg_system,cg_x,cg_b,cg_q,cg_s,cg_r,cg_r,cg_r,1e-5,0,example.number_of_CG_iterations);

        example.X+=dX;
    }


}
//#####################################################################
template class ELASTICITY_DRIVER<VECTOR<float,2> >;
template class ELASTICITY_DRIVER<VECTOR<float,3> >;
