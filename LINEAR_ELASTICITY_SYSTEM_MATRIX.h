#ifndef __LINEAR_ELASTICITY_SYSTEM_MATRIX_H__
#define __LINEAR_ELASTICITY_SYSTEM_MATRIX_H__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include "RANGE_ITERATOR.h"
#include "REGULAR_HYPERCUBE_MESH.h"

using namespace PhysBAM;

template< class T, int d >
struct LINEAR_ELASTICITY_SYSTEM_MATRIX {

    // Power Helper

    template<unsigned base, unsigned exponent> struct POWER;
    template<unsigned base> struct POWER<base,0> { enum WORKAROUND {value = 1}; };
    template<unsigned base,unsigned exponent> struct POWER { enum WORKAROUND {value=base*POWER<base,exponent-1>::value}; };

    // Types

    enum {
        VERTICES_PER_CELL = POWER<2,d>::value
    };
    typedef VECTOR<int,d> T_INDEX;
    typedef VECTOR<VECTOR<VECTOR<VECTOR<T,d>,d>,d>,d> TENSOR4; 
    typedef VECTOR<VECTOR<T,VERTICES_PER_CELL>,d> GRADIENT;
    typedef ARRAY<GRADIENT> PER_QUADRATURE_GRADIENT;
    typedef VECTOR<VECTOR<MATRIX<T,d>,VERTICES_PER_CELL>,VERTICES_PER_CELL> SYSTEM_MATRIX;
    typedef HASHTABLE<VECTOR<int,2>,MATRIX<T,d> > GLOBAL_SYSTEM_MATRIX;

    // Methods
    static void Linear_Elasticity_Tensor(const T mu, const T lambda, TENSOR4& tensor)
    {
        // 
        //  T_pkql = mu * ( e_pq * e_kl + e_pl * e_qk ) + lambda * ( e_ql * e_pk ) 
        //
        
        for(int p=1;p<=d;p++)
            for(int k=1;k<=d;k++)
                for(int q=1;q<=d;q++)
                    for(int l=1;l<=d;l++){
                        if(p==q && k==l) tensor(p)(k)(q)(l) += mu;
                        if(p==l && q==k) tensor(p)(k)(q)(l) += mu;
                        if(q==l && p==k) tensor(p)(k)(q)(l) += lambda;}
    }

    static void Gradient_Matrix( const T cell_size, VECTOR<T,d> weights, GRADIENT& G )
    {
        
        //  3D
        // 
        //  [  -(w^2)/h  -(w^2)/h      -(w^2)/h      -(w^2)/h       ((1-w)^2)/h   ((1-w)^2)/h   ((1-w)^2)/h   ((1-w)^2)/h 
        //     -(w^2)/h  -(w^2)/h       ((1-w)^2)/h   ((1-w)^2)/h  -(w^2)/h      -(w^2)/h       ((1-w)^2)/h   ((1-w)^2)/h 
        //     -(w^2)/h   ((1-w)^2)/h  -(w^2)/h       ((1-w)^2)/h  -(w^2)/h       ((1-w)^2)/h  -(w^2)/h       ((1-w)^2)/h  ]
        
        //  2D
        //
        //  [  -(w)/h  -(w)/h     (1-w)/h   (1-w)/h
        //     -(w)/h   (1-w)/h  -(w)/h     (1-w)/h  ]                     
        
        int vertex=1;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),vertex++){
            const T_INDEX& index=iterator.Index();       
            for(int v=1;v<=d;v++){
                G(v)(vertex) = 1;
                for( int w=1;w<=d;w++)
                    if(w==v)
                        if(index(w)==0)  G(v)(vertex)*=-1./cell_size;
                        else             G(v)(vertex)*=1./cell_size;
                    else
                        if(index(w)==0)  G(v)(vertex)*=(1.-weights(w));
                        else             G(v)(vertex)*=weights(w);
            }           
        }
    }
    
    static void Cell_Centered_Gradient_Matrix(const T cell_size, GRADIENT& G)
    {
        //  3D
        //
        //  [  -1/4h  -1/4h  -1/4h  -1/4h   1/4h   1/4h   1/4h   1/4h  
        //     -1/4h  -1/4h   1/4h   1/4h  -1/4h  -1/4h   1/4h   1/4h
        //     -1/4h   1/4h  -1/4h   1/4h  -1/4h   1/4h  -1/4h   1/4h   ]
        
        //  2D
        //
        //  [  -1/2h  -1/2h   1/2h   1/2h
        //     -1/2h   1/2h  -1/2h   1/2h  ]
        //   
        
        VECTOR<T,d> weights;
        weights.Fill( T(.5) );
        Gradient_Matrix( cell_size, weights, G );
    }

    static void Stiffness_Matrix(const TENSOR4& C, const PER_QUADRATURE_GRADIENT& G_q, SYSTEM_MATRIX& K)
    {
        
        //
        //  K_ij_(pq) = SUM_ C_pkql * SUM_Q ( G(Q)_ki * G(Q)_lj  ) / Q
        //
        
        T Q_max = T(G_q.m);
        
        for( int i = 1; i <= VERTICES_PER_CELL; i++ )
            for( int j = 1; j <= VERTICES_PER_CELL; j++ )
                for( int p = 1; p <= d; p++)
                    for( int q = 1; q <= d; q++)
                        for( int k = 1; k <= d; k++ )
                            for( int l = 1; l <= d; l++ ){
                                T temp = 0.0f;
                                for( int Q=1; Q<=Q_max; Q++)
                                    temp += G_q(Q)(k)(i)*G_q(Q)(l)(j);
                                K(i)(j)(p,q) += C(p)(k)(q)(l) * temp / T(Q_max);}
    }
    static void Create_Cell_Stiffness_Matrix(GLOBAL_SYSTEM_MATRIX& stiffness_matrix, const T dx, const T mu, const T lambda,const typename REGULAR_HYPERCUBE_MESH<d>::ELEMENT hypercube){
        const int number_of_quadrature_points = VERTICES_PER_CELL;   
        TENSOR4 elasticity_tensor;
        SYSTEM_MATRIX cell_stiffness_matrix;
        PER_QUADRATURE_GRADIENT elasticity_gradient_matrices;
        VECTOR<VECTOR<T,d>,number_of_quadrature_points> weights;
        
        elasticity_gradient_matrices.Resize(number_of_quadrature_points);
        int quadrature_point=1;
        for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(T_INDEX(),T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next(),quadrature_point++){
            const T_INDEX& index=iterator.Index();
            for(int i=1;i<=d;i++)
                if(index(i)==0)  weights(quadrature_point)(i)=(1-one_over_root_three)*.5;
                else             weights(quadrature_point)(i)=(1+one_over_root_three)*.5;
            Gradient_Matrix( dx, weights(quadrature_point), elasticity_gradient_matrices(quadrature_point) );}

        Linear_Elasticity_Tensor(mu,lambda,elasticity_tensor);
        Stiffness_Matrix(elasticity_tensor, elasticity_gradient_matrices,cell_stiffness_matrix);

        int vertex1=1;
        for(RANGE_ITERATOR<d> node_iterator1(RANGE<T_INDEX>::Unit_Box());node_iterator1.Valid();node_iterator1.Next(),++vertex1){
            //const T_INDEX& index1=node_iterator1.Index();
            const int node1=hypercube(vertex1);
            int vertex2=1;
            for(RANGE_ITERATOR<d> node_iterator2(RANGE<T_INDEX>::Unit_Box());node_iterator2.Valid();node_iterator2.Next(),++vertex2){
                //const T_INDEX& index2=node_iterator2.Index();
                const int node2=hypercube(vertex2);
                // Is the cell stiffness matrix symetric?
                // errr.... should it be minus equal here? Otherwise the CG will complaint....
                stiffness_matrix.Get_Or_Insert(VECTOR<int,2>(node1,node2),MATRIX<T,d>())-=cell_stiffness_matrix(vertex1)(vertex2)*dx*dx*dx;}}
    }
};



#endif
