//#####################################################################
// Copyright 2014, Nathan Mitchell
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REGULAR_HYPERCUBE_MESH
//#####################################################################
#ifndef __REGULAR_HYPERCUBE_MESH__
#define __REGULAR_HYPERCUBE_MESH__

// #include <PhysBAM_Tools/Log/LOG.h>
// #include <PhysBAM_Tools/Arrays/ARRAY.h>
// #include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
// #include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
// #include <PhysBAM_Tools/Data_Structures/PAIR.h>
// #include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
// #include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
// #include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>

namespace PhysBAM{ 

// namespace POWER{
// template<unsigned base, unsigned exponent> struct POWER;
// template<unsigned base> struct POWER<base,0> { enum WORKAROUND {value = 1}; };
// template<unsigned base,unsigned exponent> struct POWER { enum WORKAROUND {value=base*POWER<base,exponent-1>::value}; };
// }


template<int d>
class REGULAR_HYPERCUBE_MESH
{
    
    typedef VECTOR<int,d> T_INDEX;

public:
    enum { vertices_per_cell=1<<d };
    enum { max_neighbors=d*2 };
    typedef VECTOR<int,vertices_per_cell> ELEMENT;

 
    ARRAY<ELEMENT> elements;
    ARRAY<T_INDEX> element_coordinates;
    ARRAY<ARRAY<int> > *neighbors;
    
 
    REGULAR_HYPERCUBE_MESH():neighbors(NULL) {}
    ~REGULAR_HYPERCUBE_MESH() {
         if( neighbors )
             delete neighbors;
         neighbors = NULL;
     }

    void Compute_Neighbors(){
         if( neighbors )
             delete neighbors;
         neighbors = new ARRAY<ARRAY<int> >(elements.m);

         HASHTABLE<T_INDEX,int > cell_hash;
         for (int i=1;i<=elements.m;i++) {
            cell_hash.Insert(element_coordinates(i),i);
         }

         VECTOR<VECTOR<int,d>, 2> offset;
         offset(1)=T_INDEX::All_Ones_Vector();
         offset(2)=T_INDEX::All_Ones_Vector()*(-1);
         for (int i=1;i<=elements.m;i++) {
             T_INDEX cell_index=element_coordinates(i);
             for (int v=1;v<=d;v++) {
                for (int u=1;u<=2;u++) {
                    T_INDEX potential_neighbor=cell_index;
                    potential_neighbor(v)+=offset(u)(v);
                    if (cell_hash.Contains(potential_neighbor))
                        (*neighbors)(i).Append(cell_hash.Get(potential_neighbor));
                }
             }
         }

    //        vertices_per_face=1<<(d-1),
    //        faces_per_cell=2*d
    // };

    // typedef VECTOR<T,d> TV;
    // typedef VECTOR<int,d> T_INDEX;
    // typedef VECTOR<int,vertices_per_face> T_FACE;


    // struct T_CELL {
    //     VECTOR<int, vertices_per_cell> vertices;
    //     T_INDEX index;            

    //     static T_FACE face_indices(int axis, bool direction) {
    //         VECTOR<int, vertices_per_face> f;
    //         int flat_index=1;
    //         int face_index=1;
    //         for(RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1)); iter.Valid(); iter.Next()){
    //             const T_INDEX& index=iter.Index();
    //             if( index(axis) == direction ){
    //                 f(face_index) = flat_index;
    //                 face_index++;
    //             }
    //             flat_index++;
    //         }
    //         PHYSBAM_ASSERT(face_index == vertices_per_face+1);
    //         return f;
    //     };

    //     T_FACE face(int axis, bool direction) const {
    //         VECTOR<int, vertices_per_face> f;
    //         int flat_index=1;
    //         int face_index=1;
    //         for(RANGE_ITERATOR<d> iter(RANGE<T_INDEX>(T_INDEX(),T_INDEX()+1)); iter.Valid(); iter.Next()){
    //             const T_INDEX& index=iter.Index();
    //             if( index(axis) == direction ){
    //                 f(face_index) = vertices(flat_index);
    //                 face_index++;
    //             }
    //             flat_index++;
    //         }
    //         PHYSBAM_ASSERT(face_index == vertices_per_face+1);
    //         return f;
    //     };

    //     VECTOR<T_FACE,faces_per_cell> faces() const {
    //         VECTOR<T_FACE,faces_per_cell> f;
    //         for( int i=0; i<d; i++)
    //             for( int a=0; a<2; a++)
    //                 f( (2*i + a) +1 ) = face(i+1, a);
    //         return f;
    //     };    
           
    // };

    // Members
 
    // T dx;
    // RANGE<TV> domain;

    // // Acceleration Structures

    // // Constructors
        
    // REGULAR_HYPERCUBE_MESH() : dx(0.0), neighbors(NULL) {};
        
    // REGULAR_HYPERCUBE_MESH(const T dx_input) : dx(dx_input), neighbors(NULL) {};
        
    // REGULAR_HYPERCUBE_MESH(const REGULAR_HYPERCUBE_MESH& mesh) : dx(mesh.dx), elements(mesh.elements), 
    //                                                             neighbors(NULL) {};

    // template<class GRID>
    // REGULAR_HYPERCUBE_MESH(const GRID& grid) : dx(grid.dX.Max()), domain(grid.domain), neighbors(NULL)
    // {
    //     HASHTABLE<T_INDEX, int> node_map;
    //     int max_node_plus_one = 1;

    //     for(RANGE_ITERATOR<d> cell_iter(grid.Cell_Indices()); cell_iter.Valid(); cell_iter.Next()){
    //         T_CELL new_cell;
    //         new_cell.index = cell_iter.Index();
    //         int vpos = 1;
    //         for(RANGE_ITERATOR<d> node_iter(RANGE<T_INDEX>(cell_iter.Index(),cell_iter.Index()+1)); node_iter.Valid(); node_iter.Next(), vpos++){
    //             int vertex = node_map.Get_Or_Insert( node_iter.Index(), max_node_plus_one );
    //             if(vertex == max_node_plus_one)
    //                 max_node_plus_one++;
    //             new_cell.vertices(vpos) = vertex;                    
    //         }
    //         elements.Append(new_cell);
    //     }
            
    // };

    // REGULAR_HYPERCUBE_MESH( const T dx_input, const RANGE<TV>& box ):
    //     dx(dx_input), neighbors(NULL)
    // {
    //     T_INDEX cell_counts = T_INDEX( ceil(box.Edge_Lengths()/dx) );
    //     HASHTABLE<T_INDEX, int> node_map;
    //     int max_node_plus_one = 1;

    //     RANGE<T_INDEX> cell_range( T_INDEX()+1, cell_counts );

    //     domain.min_corner = box.min_corner;
    //     domain.max_corner = box.min_corner + TV(cell_range.max_corner)*dx;

    //     for(RANGE_ITERATOR<d> cell_iter(cell_range); cell_iter.Valid(); cell_iter.Next()){
    //         T_CELL new_cell;
    //         new_cell.index = cell_iter.Index();
    //         int vpos = 1;
    //         for(RANGE_ITERATOR<d> node_iter(RANGE<T_INDEX>(cell_iter.Index(),cell_iter.Index()+1)); node_iter.Valid(); node_iter.Next(), vpos++){
    //             int vertex = node_map.Get_Or_Insert( node_iter.Index(), max_node_plus_one );
    //             if(vertex == max_node_plus_one)
    //                 max_node_plus_one++;
    //             new_cell.vertices(vpos) = vertex;                    
    //         }
    //         elements.Append(new_cell);
    //     } 
            
    // };

 
        

    // // Methods
    // TV Node( const int cell, const int v ) const{
    //     T_INDEX offset;
    //     for( int i=0; i<d; i++)
    //         offset(d-i) = int( ((v-1) & 1<<i) > 0 );
    //     return Node(elements(cell).index + offset);
    // }

    // TV Node( const T_INDEX& index) const{
    //     return domain.min_corner + TV(index-1) * dx;
    // }

    // TV Cell_Center( const int& cell) const{
    //     const T_INDEX& index = elements(cell).index;
    //     return Cell_Center(index);
    // }

    // TV Cell_Center( const T_INDEX& index) const{
    //     return domain.min_corner+(TV(index)+(T).5)*dx;
    // }
        
    // // Acceleration Structure Initializers
      
 








        //   HASHTABLE<T_FACE, int> face_map;
        //   HASHTABLE<int, TRIPLE<int,int,int> > face_ownership;
        //   int max_faces_plus_one = 1;
        //   for( int c = 1; c <= elements.m; c ++ ){
        //       VECTOR<T_FACE,faces_per_cell> fs = elements(c).faces();
        //       for( int i=1; i <= faces_per_cell; i++){
        //           int t = face_map.Get_Or_Insert(fs(i),max_faces_plus_one);
        //           if(t==max_faces_plus_one)
        //               max_faces_plus_one++;
        //           TRIPLE<int,int,int>& owner_cells = face_ownership.Get_Or_Insert( t, TRIPLE<int,int,int>((i-1)/2,-1,-1) );
        //           PHYSBAM_ASSERT( owner_cells.x == (i-1)/2 );
        //           if( (i-1)%2 ){
        //               PHYSBAM_ASSERT( owner_cells.y == -1 || owner_cells.y == c );
        //               owner_cells.y = c;
        //           }
        //           else{
        //               PHYSBAM_ASSERT( owner_cells.z == -1 || owner_cells.z == c );
        //               owner_cells.z = c;                                            
        //           }
        //       }
        //   }

        //   neighbors->Resize( elements.m );
        //   for( HASHTABLE_ITERATOR<int, TRIPLE<int,int,int> > iterator( face_ownership ); iterator.Valid(); iterator.Next() ) {
        //    
        //       TRIPLE<int,int,int> data = iterator.Data();
        //       if(data.y != -1) (*neighbors)(data.y)(data.x+1).y = data.z;
        //       if(data.z != -1) (*neighbors)(data.z)(data.x+1).x = data.y;
        //   }
        //    
     }


};

}

#endif
