//#####################################################################
// Copyright 2014-2016, Michael Doescher, Haixiang Liu, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_SURFACE
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "EMBEDDED_SURFACE.h"
#include "LINEAR_ELASTICITY_SYSTEM_MATRIX.h"
#include "EMBEDDING_TOOLS.h"

#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>

#include "RANGE_ITERATOR.h"
#include "REGULAR_HYPERCUBE_MESH.h"
#include "MATERIAL_MODEL.h"
#include "COROTATED.h"

using namespace PhysBAM;

struct COROTATED_TAG;

//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Initialize(){
    BASE::Initialize();
   
    // example specific model + domain resizing 
    Load_Embedded_Surface();
    TV min_corner=domain.Minimum_Corner();

    // initialize mesh
    for(RANGE_ITERATOR<d> cell_iterator(RANGE<T_INDEX>(T_INDEX(),size-1));cell_iterator.Valid();cell_iterator.Next()){
        const T_INDEX& cell_index=cell_iterator.Index();
        // Normal element insertion
        T_ELEMENT element;
        int vertex=1;
        for(RANGE_ITERATOR<d> node_iterator(RANGE<T_INDEX>(cell_index,cell_index+1));node_iterator.Valid();node_iterator.Next(),vertex++){
            const T_INDEX& node_index=node_iterator.Index();
            element(vertex)=node_hash.Get_Or_Insert(node_index,node_hash.Size()+1);
            node_to_ijk.Get_Or_Insert( element(vertex), node_index );}
        mesh.elements.Append(element);
        mesh.element_coordinates.Append(cell_index);
    }
    mesh.Compute_Neighbors();
       
    // initialize X 
    X.Resize(node_hash.Size());
    for(HASHTABLE_ITERATOR<T_INDEX,const int> iterator(node_hash);iterator.Valid();iterator.Next())
            X(iterator.Data())=(TV(iterator.Key())*dx) + min_corner;
    X_resting = X;

    // rasterize and map embedded vertices then select example specific kinematic nodes
    EMBEDDINGTOOLS<T,d>::Rasterize(embedded_vertices, embedded_faces, dx, mesh, X, cell_is_full);   
    EMBEDDINGTOOLS<T,d>::Generate_Embedding_Map(embedded_vertices, dx, mesh, X, embedding_map);   
    Select_Kinematic_Nodes();
}




//#####################################################################
// Function Clear_Values_Of_Kinematic_Particles
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const{
    for(HASHTABLE_ITERATOR<int> iterator(node_is_dirichlet); iterator.Valid(); iterator.Next()){
        array(iterator.Key())=TV();
    }   
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Register_Options(){

    BASE::Register_Options();
    if (d==2) {
        parse_args->Add_Vector_2D_Argument("-size",VECTOR<double,2>(80,72),"n, n","Size of the cell grid");
        parse_args->Add_Vector_2D_Argument("-macroblock_size",VECTOR<double,2>(16,8),"b b","Size of the macro blocks");
    }
    if (d==3) {
        parse_args->Add_Vector_3D_Argument("-size",VECTOR<double,3>(80,72,72),"n n n","Size of the cell grid");
        parse_args->Add_Vector_3D_Argument("-macroblock_size",VECTOR<double,3>(16,8,8),"b b b","Size of the macro blocks");
    }
    parse_args->Add_Integer_Argument("-example",3,"n","Which test example");
    parse_args->Add_Double_Argument("-rate",0.05,"d","Rate of kinematic object movement in m/s (note: frame_rate independent)");
    parse_args->Add_Double_Argument("-perturbation_magnitude",0.1,"d","for randomly perturbing node positions");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Parse_Options(){
    BASE::Parse_Options();
    for(int axis=1;axis<=d;++axis){ 
        if(d==2) size(axis)=parse_args->Get_Vector_2D_Value("-size")(axis);
        if(d==3) size(axis)=parse_args->Get_Vector_3D_Value("-size")(axis);}
    for(int axis=1;axis<=d;++axis) {
        if(d==2) macroblock_size(axis)=parse_args->Get_Vector_2D_Value("-macroblock_size")(axis);
        if(d==3) macroblock_size(axis)=parse_args->Get_Vector_3D_Value("-macroblock_size")(axis);
    }
    example=parse_args->Get_Integer_Value("-example");
    rate=parse_args->Get_Double_Value("-rate");
    rate/=frame_rate;
    perturbation_magnitude=parse_args->Get_Double_Value("-perturbation_magnitude");
}
//#####################################################################
// Function Add_Elastic_Forces
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const
{
    typedef LINEAR_ELASTICITY_SYSTEM_MATRIX<T,d> LESM;
    typedef VECTOR<int,d> T_INDEX;

    VECTOR< TV,  REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell> quadrature_weights;
    {int qp=1;
    for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));
        iterator.Valid();iterator.Next(),qp++){
        const T_INDEX& index=iterator.Index();
        for(int i=1;i<=d;i++)
            if(index(i)==0)  quadrature_weights(qp)(i)=(1-one_over_root_three)*.5;
            else             quadrature_weights(qp)(i)=(1+one_over_root_three)*.5;
    }}
        
    for( int element=1; element <= mesh.elements.m; element++){
        if (!(cell_is_full.Contains(element))) continue;

        MATRIX_MXN<T> Ds;Ds.Resize(d,REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell);
        MATRIX_MXN<T> H;H.Resize(d,REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell);
        T Ds_c[3][8];
        T H_c[3][8];
        {
            for(int n =1; n <= REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell; n++ ){
                const int node = mesh.elements(element)(n);
                for(int v=1;v<=d;v++){
                    Ds(v,n)=X(node)(v);
                    Ds_c[v-1][n-1]=X(node)(v);
                    H_c[v-1][n-1] = 0.0f;
                }
            }
        }
        
        
// #define USE_KERNEL_GRADIENTS
// #define USE_KERNEL_SVD
// #define USE_KERNEL_P

        for( int qp=1; qp <= REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell; qp ++ ){ 
            MATRIX<T,d> Fe;

#if !defined(USE_KERNEL_GRADIENTS)
            typename LESM::GRADIENT G_tmp;
            LESM::Gradient_Matrix( dx, quadrature_weights(qp), G_tmp);                     

            MATRIX_MXN<T> G; G.Resize(d,REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell);
            for(int n =1; n <= REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell; n++ )
                for(int v=1;v<=d;v++){
                    G(v,n) = G_tmp(v)(n);                    
                }

            // Don't add one, uses exact positions already
            Fe=MATRIX<T,d>(Ds*G.Transposed());
#else
            T W_c[d];
            for( int w=1; w<=d; w++)
                W_c[w-1] = quadrature_weights(qp)(w);
            T one_over_h = 1.0/dx;
            T F_c[9];
            Weighted_Gradient<float,float,int>( Ds_c, F_c, W_c, one_over_h );            
            for( int x=0; x<9;x++)
                Fe.x[x]=F_c[x];
#endif


            MATRIX<T,d> Ue,Ve;
            DIAGONAL_MATRIX<T,d> Sigmae;
#if !defined(USE_KERNEL_SVD)
            Fe.Fast_Singular_Value_Decomposition(Ue,Sigmae,Ve);
#else
            T U_c[9], V_c[9];
            T Sigma_c[3];
            Singular_Value_Decomposition<float,float,int>( F_c, U_c, Sigma_c, V_c );
            for( int i=0; i<9; i++){
                Ue.x[i] = U_c[i];
                Ve.x[i] = V_c[i];
                if(i<3)
                    Sigmae(i+1) = Sigma_c[i];
            }            
#endif
            Sigmae=Sigmae.Clamp_Min(1e-4);

            MATRIX<T,d> P;
#if !defined(USE_KERNEL_P)
            DIAGONAL_MATRIX<T,d> P_hat;
            // Corotated P Formula
            P_hat = 2*constant_mu*(Sigmae-1) + constant_lambda*(Sigmae-1).Trace();
            P=Ue*P_hat.Times_Transpose(Ve);
#else
            Add_Force_Single_QPoint<COROTATED_TAG,float,float,int>::Run( constant_mu, 0, 0, constant_lambda, U_c, V_c, Sigma_c, F_c );
            for( int i=0; i<9; i++)
                P.x[i] = F_c[i]; // Copy out final P
#endif

            T cell_volume;
            if(d==2)
                cell_volume = dx*dx;
            else
                cell_volume = dx*dx*dx;

#if !defined(USE_KERNEL_GRADIENTS)
            H += -P*G*cell_volume*(1.0 / REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell);
#else          
            for( int x=0; x<9;x++)
                F_c[x]=P.x[x];
            Weighted_Accumulation<float,float,int>( H_c, F_c, W_c, one_over_h, (-1.0f * dx*dx*dx / REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell) ); 
#endif     

        }

#if defined(USE_KERNEL_GRADIENTS)
        for(int n =1; n <= REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell; n++ )
            for(int v=1;v<=d;v++){
                H(v,n) = H_c[v-1][n-1];
            }
#endif

        {
            for(int n =1; n <= REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell; n++ ){
                const int node = mesh.elements(element)(n);
                for(int v=1;v<=d;v++)
                    force(node)(v) += H(v,n);
            }
        }

        
    }

}
//==========================================================
//            BUILD  M  Matrix
//==========================================================

template<class T>
MATRIX<T,2> Build_M(const ROTATED_STRESS_DERIVATIVE<T,2>& dPdF,const VECTOR<T,2>& hi,const VECTOR<T,2>& hj)
{
    MATRIX<T,2> M;

//     M(1,1)=dPdF.a1111*hi(1)*hj(1)
//           +dPdF.a1212*hi(2)*hj(2)
//           +dPdF.a1313*hi(3)*hj(3);
// 
//     M(2,1)=dPdF.a1122*hi(2)*hj(1)
//           +dPdF.a1221*hi(1)*hj(2);
// 
//     M(3,1)=dPdF.a1133*hi(3)*hj(1)
//           +dPdF.a1331*hi(1)*hj(3);
// 
//     M(1,2)=dPdF.a1122*hi(1)*hj(2)
//           +dPdF.a1221*hi(2)*hj(1);
// 
//     M(2,2)=dPdF.a1212*hi(1)*hj(1)
//           +dPdF.a2222*hi(2)*hj(2)
//           +dPdF.a2323*hi(3)*hj(3);
// 
//     M(3,2)=dPdF.a2233*hi(3)*hj(2)
//           +dPdF.a2332*hi(2)*hj(3);
// 
//     M(1,3)=dPdF.a1133*hi(1)*hj(3)
//           +dPdF.a1331*hi(3)*hj(1);
// 
//     M(2,3)=dPdF.a2233*hi(2)*hj(3)
//           +dPdF.a2332*hi(3)*hj(2);
// 
//     M(3,3)=dPdF.a1313*hi(1)*hj(1)
//           +dPdF.a2323*hi(2)*hj(2)
//           +dPdF.a3333*hi(3)*hj(3);
// 
    return M;
}
template<class T>
MATRIX<T,3> Build_M(const ROTATED_STRESS_DERIVATIVE<T,3>& dPdF,const VECTOR<T,3>& hi,const VECTOR<T,3>& hj)
{
    MATRIX<T,3> M;

    M(1,1)=dPdF.a1111*hi(1)*hj(1)
          +dPdF.a1212*hi(2)*hj(2)
          +dPdF.a1313*hi(3)*hj(3);

    M(2,1)=dPdF.a1122*hi(2)*hj(1)
          +dPdF.a1221*hi(1)*hj(2);

    M(3,1)=dPdF.a1133*hi(3)*hj(1)
          +dPdF.a1331*hi(1)*hj(3);

    M(1,2)=dPdF.a1122*hi(1)*hj(2)
          +dPdF.a1221*hi(2)*hj(1);

    M(2,2)=dPdF.a1212*hi(1)*hj(1)
          +dPdF.a2222*hi(2)*hj(2)
          +dPdF.a2323*hi(3)*hj(3);

    M(3,2)=dPdF.a2233*hi(3)*hj(2)
          +dPdF.a2332*hi(2)*hj(3);

    M(1,3)=dPdF.a1133*hi(1)*hj(3)
          +dPdF.a1331*hi(3)*hj(1);

    M(2,3)=dPdF.a2233*hi(2)*hj(3)
          +dPdF.a2332*hi(3)*hj(2);

    M(3,3)=dPdF.a1313*hi(1)*hj(1)
          +dPdF.a2323*hi(2)*hj(2)
          +dPdF.a3333*hi(3)*hj(3);

    return M;
}

//#####################################################################
// Function Update Position Based State
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Update_Position_Based_State(T time){
    typedef LINEAR_ELASTICITY_SYSTEM_MATRIX<T,d> LESM;
    typedef VECTOR<int,d> T_INDEX;

    stiffness_matrix.Clean_Memory();
    
    enum {vertices_per_cell=1<<d, number_of_quadrature_points=1<<d};

    VECTOR< TV,  number_of_quadrature_points> quadrature_weights;
    {int qp=1;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));
            iterator.Valid();iterator.Next(),qp++){
            const T_INDEX& index=iterator.Index();
            for(int i=1;i<=d;i++)
                if(index(i)==0)  quadrature_weights(qp)(i)=(1-one_over_root_three)*.5;
                else             quadrature_weights(qp)(i)=(1+one_over_root_three)*.5;
        }}
    
    typename LESM::GRADIENT G_tmp;       
    VECTOR< MATRIX_MXN<T>, number_of_quadrature_points > G;
    for( int i=1; i<=number_of_quadrature_points; i++) 
        G(i).Resize(d,REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell);
    
    for( int qp=1; qp <= number_of_quadrature_points; qp ++ ){ 
        LESM::Gradient_Matrix( dx, quadrature_weights(qp), G_tmp);                     
        for(int n =1; n <= REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell; n++ )
            for(int v=1;v<=d;v++){
                G(qp)(v,n) = G_tmp(v)(n);                    
            }
    }
    

    for(int element=1;element<=mesh.elements.m;element++){
        if (!(cell_is_full.Contains(element))) continue;
        MATRIX_MXN<T> Ds;Ds.Resize(d,REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell);
        {
            for(int n =1; n <= REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell; n++ ){
                const int node = mesh.elements(element)(n);
                for(int v=1;v<=d;v++){
                    Ds(v,n)=X(node)(v);
                }
            }
        }

        for( int qp=1; qp <= number_of_quadrature_points; qp ++ ){ 
            MATRIX<T,d> Fe;
            Fe=Ds.Times_Transpose(G(qp));
            
            MATRIX<T,d> Ue,Ve;
            DIAGONAL_MATRIX<T,d> Sigmae;
            Fe.Fast_Singular_Value_Decomposition(Ue,Sigmae,Ve);
            Sigmae=Sigmae.Clamp_Min(1e-4);
            
            // This is only identity for corotated elasticity
            ROTATED_STRESS_DERIVATIVE<T,d> dP_dFe;
            dP_dFe=MATERIAL_MODEL<COROTATED<T,d> >::Rotated_Stress_Derivative(Sigmae,constant_mu,
                                                                              std::min((T)5.*constant_mu,constant_lambda));
            dP_dFe.Make_Positive_Definite();

            MATRIX_MXN<T> He;
            He.Resize(d,REGULAR_HYPERCUBE_MESH<d>::vertices_per_cell);
            He = Ve.Transpose_Times(G(qp));

            T cell_volume;
            if(d==2)
                cell_volume = dx*dx;
            else
                cell_volume = dx*dx*dx;
            T W0 = -cell_volume * (1.0f / number_of_quadrature_points);

            for(int i=1; i<=vertices_per_cell; i++)
                for(int j=1; j<=vertices_per_cell; j++){
                    
                    VECTOR<T,d> hi, hj;
                    for( int q=1; q<=d; q++){
                        hi(q) = He(q,i);
                        hj(q) = He(q,j);
                    }
                    MATRIX<T,d> M = Build_M( dP_dFe, hi, hj );
                    M = M*W0;                    
                    stiffness_matrix.Get_Or_Insert(VECTOR<int,2>(mesh.elements(element)(i), mesh.elements(element)(j)), MATRIX<T,d>()) += (Ue)*(M)*(Ue.Transposed());
                }



        }
    }
}

//#####################################################################
// template class EMBEDDED_SURFACE<VECTOR<float,2> >;
 template class EMBEDDED_SURFACE<VECTOR<float,3> >;
