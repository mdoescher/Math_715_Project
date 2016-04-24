#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>

#include "EMBEDDING_TOOLS.h"


using namespace PhysBAM;

//#####################################################################
// Function Get_Grid_Bounds
//#####################################################################
template< class T, int d>
void EMBEDDINGTOOLS<T,d>::Get_Grid_Bounds( const ARRAY<VECTOR<T,d> >& vertices, const ARRAY<VECTOR<int,d> >& triangles, const T dx, VECTOR<int, d>& cell_bounds, VECTOR<T, d>& min_corner)
{
    LOG::SCOPE scope("Embedding_Tools::Get_Grid_Bounds()");

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    
    // Compute (padded) bounding box
    
    RANGE<TV> box=RANGE<TV>::Empty_Box();
    for(int v=1;v<=vertices.m;v++)
                box.Enlarge_To_Include_Point(vertices(v));
    box.Change_Size(.05*dx);
    
    // Subdivide box (potentially enlarging it) so that the cells along each dimension == 3 (mod 4)
    
    for(int v=1;v<=3;v++)
        cell_bounds(v)=((1+(int)(box.Edge_Lengths()(v)/dx))/4)*4+3;
    box.max_corner=box.min_corner+TV(cell_bounds)*dx;        
    min_corner = box.min_corner;              
    return;
}


//#####################################################################
// Function Generate_Embedding_Map (macroblock version)
//#####################################################################
template<class T, int d> void EMBEDDINGTOOLS<T,d>::
Generate_Embedding_Map(const ARRAY<VECTOR<T,d> > &vertices, const T dx, 
        const REGULAR_HYPERCUBE_MESH<d>& mesh, const ARRAY<VECTOR<T,d> >& X, ARRAY<PAIR<int,VECTOR<T,d> > >& embedding_map)
{
    // note: ARRAY<PAIR<int,VECTOR<T,d> > >& embedding_map for each embedded node PAIR(cell, multilinear coordinates)
    LOG::SCOPE scope("Embedding_Tools::Generate_Embedding_Map()");
    embedding_map.Resize(vertices.m);

    typedef VECTOR<T,d> TV;
    // generate box hierarcy for the MESH    
    ARRAY<RANGE<TV> > boxes;
    for (int i=1;i<=mesh.elements.m;i++) {
        RANGE<TV> box=RANGE<TV>::Empty_Box();
        for (int j=1;j<=mesh.vertices_per_cell;j++) 
            box.Enlarge_To_Include_Point(X(mesh.elements(i)(j)));
        boxes.Append(box);}
    BOX_HIERARCHY<TV> box_hierarchy;
    box_hierarchy.Set_Leaf_Boxes(boxes,true);

    // traverse the list of embedded_vertices and identify which cells are intersected by each vertex
    for(int i=1;i<=vertices.m;i++){
        RANGE<TV> vertex_box=RANGE<TV>::Empty_Box();
        vertex_box.Enlarge_To_Include_Point(vertices(i));
        ARRAY<int> intersection_list;
        box_hierarchy.Intersection_List(vertex_box, intersection_list, (T)1e-5);
        TV coordinates;
        int min_node, max_node;
        int cell=0;
        for (int j=1;j<=intersection_list.m;j++) {
            cell=intersection_list(j);
            min_node=mesh.elements(cell)(1), max_node=mesh.elements(cell)(8);
//            for (int j=2;j<=mesh.vertices_per_cell;j++) {
//                if( (X(mesh.elements(cell)(j))-X(min_node)).Min() < (T)-1e-4 ) min_node=mesh.elements(cell)(j);
//                if( (X(max_node)-X(mesh.elements(cell)(j))).Min() < (T)-1e-4 ) max_node=mesh.elements(cell)(j);
//            }
            coordinates=(vertices(i)-X(min_node))/(X(max_node)-X(min_node));
            if (coordinates.Min() >= T()) break;
        }
        TV clamped_coordinates=clamp_min(coordinates,TV());
        embedding_map(i)=PAIR<int,TV>(cell,clamped_coordinates);
    }
}


//#####################################################################
// Function Rasterize (macroblock version)
//#####################################################################
template<class T, int d> void EMBEDDINGTOOLS<T,d>::
Rasterize(const ARRAY<VECTOR<T,d> > &vertices, const ARRAY<VECTOR<int,d> > &faces, const T dx, 
        const REGULAR_HYPERCUBE_MESH<d>& mesh, const ARRAY<VECTOR<T,d> >& X, HASHTABLE<int>& cell_is_full) 
{
 LOG::SCOPE scope("Embedding_Tools::Rasterize()");

    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> T_INDEX;
    
    {
        LOG::SCOPE scope("Marking all grid cells intersecting mesh.");

        // generate box hierarcy for the MESH    
        ARRAY<RANGE<TV> > boxes;
        for (int i=1;i<=mesh.elements.m;i++) {
            RANGE<TV> box=RANGE<TV>::Empty_Box();
            for (int j=1;j<=mesh.vertices_per_cell;j++) 
                box.Enlarge_To_Include_Point(X(mesh.elements(i)(j)));
            boxes.Append(box);}
        BOX_HIERARCHY<TV> box_hierarchy;
        box_hierarchy.Set_Leaf_Boxes(boxes,true);

        // traverse the list of triangles and identify which cells are intersected by each triangle
        for(int t=1;t<=faces.m;t++){
            RANGE<TV> triangle_box=RANGE<TV>::Empty_Box();
            for (int v=1;v<=d;v++) triangle_box.Enlarge_To_Include_Point(vertices(faces(t)(v)));
            ARRAY<int> intersection_list;
            box_hierarchy.Intersection_List(triangle_box, intersection_list);
            TRIANGLE_3D<T> triangle(vertices(faces(t)(1)), vertices(faces(t)(2)), vertices(faces(t)(3)));
            for (int i=1;i<=intersection_list.m;i++) {
                RANGE<TV> cell = boxes(intersection_list(i));   
                if(INTERSECTION::Intersects<T>(cell,triangle,(T)1e-4)) cell_is_full.Set(intersection_list(i));
            }
        }
    }
    
     // Flood fill the exterior region up to the boundary
     HASHTABLE<int> cell_is_exterior;
     {
         LOG::SCOPE scope("Flood filling all nodes to determine outside region.");
         ARRAY<bool> cell_added_to_stack(mesh.elements.m);
         STACK<int> stack;
         // add all domain boundary nodes to the stack (they will either be already filled by the embedded model or exterior)
         for (int i=1;i<=mesh.elements.m;i++) 
             if ((*mesh.neighbors)(i).m < mesh.max_neighbors) {
                 stack.Push(i);
                 cell_added_to_stack(i)=true;
             }

         while(!stack.Empty()){
             const int cell=stack.Pop();
             if (cell_is_full.Contains(cell)) continue;
             cell_is_exterior.Set(cell);
             for (int i=1;i<=(*mesh.neighbors)(cell).m;i++) {
                 int neighbor_cell=(*mesh.neighbors)(cell)(i);
                 if (cell_added_to_stack(neighbor_cell)==false) {
                     stack.Push(neighbor_cell);
                     cell_added_to_stack(neighbor_cell)=true;
                 }
             }
         }
     }

     // Mark remaining cells as full
     {
         LOG::SCOPE scope("Marking all remaining cells as inside.");
         {
             for (int i=1;i<=mesh.elements.m;i++) 
                 if (!cell_is_exterior.Contains(i)) cell_is_full.Set(i);
         }
     }
}

//#####################################################################
// Function Coarsen
//#####################################################################

    template<class T, int d>
void EMBEDDINGTOOLS<T,d>::Coarsen( const GRID<VECTOR<T,d> >& fine_domain, const RANGE<VECTOR<int,d> >& fine_unpadded_domain,
        const RANGE<VECTOR<int,d> >& fine_padded_domain, const ARRAY< bool, VECTOR<int, d> >& fine_voxmap,
        const GRID<VECTOR<T,d> >& coarse_domain, const RANGE<VECTOR<int,d> >& coarse_unpadded_domain,
        const RANGE<VECTOR<int,d> >& coarse_padded_domain, ARRAY< CELL_TYPE , VECTOR<int, d> >& coarse_voxmap)
 {
     LOG::SCOPE scope("Embedding_Tools::Coarsen()");
 
     typedef VECTOR<T,d> TV;
     typedef VECTOR<int,d> T_INDEX;
     
     int interior_count=0;
     int exterior_count=0;
     int boundary_count=0;
 
     // We start at the same place.
     PHYSBAM_ASSERT( fine_domain.Xmin() == coarse_domain.Xmin() );
     int coarsen_ratio = coarse_domain.Maximum_Edge_Length() / fine_domain.Maximum_Edge_Length();
     
     coarse_voxmap.Fill( EXTERIOR_CELL_TYPE );
     for( RANGE_ITERATOR<d> coarse_iterator( coarse_unpadded_domain ); coarse_iterator.Valid(); coarse_iterator.Next() ){
         bool is_full = true;
         bool is_empty = true;
         
         const T_INDEX& index = coarse_iterator.Index();
         T_INDEX fine_start = (index-1)*coarsen_ratio+1;
         T_INDEX fine_end = (index-1)*coarsen_ratio+T_INDEX::All_Ones_Vector()*coarsen_ratio;
 
         for( RANGE_ITERATOR<d> fine_iterator( RANGE<T_INDEX>(fine_start, fine_end)); fine_iterator.Valid(); fine_iterator.Next()){
                 if( fine_voxmap( fine_iterator.Index() ) == true)
                     is_empty = false;
                 if( fine_voxmap( fine_iterator.Index() ) == false)
                     is_full = false;
             }
             
         if( is_empty ){
             coarse_voxmap( index ) = EXTERIOR_CELL_TYPE;
             exterior_count++;
         }
         if( is_full ){
             coarse_voxmap( index ) = INTERIOR_CELL_TYPE;
             interior_count++;
         }
         if( !is_full && !is_empty ){
             coarse_voxmap( index ) = INTERIOR_CELL_TYPE; 
             interior_count++;
         }
     }
 
     LOG::cout << "Coarse Grid has " << interior_count << " INTERIOR_CELLS" << std::endl;
     LOG::cout << "Coarse Grid has " << exterior_count << " EXTERIOR_CELLS" << std::endl;
     LOG::cout << "Coarse Grid has " << boundary_count << " BOUNDARY_CELLS" << std::endl;
 }
 
 
 //#####################################################################
 // Function Coarsen
 //#####################################################################
 
 template<class T, int d>
 void EMBEDDINGTOOLS<T,d>::Coarsen_Density( const GRID<VECTOR<T,d> >& fine_domain, const RANGE<VECTOR<int,d> >& fine_unpadded_domain,
                               const RANGE<VECTOR<int,d> >& fine_padded_domain, const ARRAY< bool, VECTOR<int, d> >& fine_voxmap,
                               const GRID<VECTOR<T,d> >& coarse_domain, const RANGE<VECTOR<int,d> >& coarse_unpadded_domain,
                               const RANGE<VECTOR<int,d> >& coarse_padded_domain, ARRAY< float , VECTOR<int, d> >& coarse_voxmap)
 {
     LOG::SCOPE scope("Embedding_Tools::Coarsen_Density()");
 
     typedef VECTOR<T,d> TV;
     typedef VECTOR<int,d> T_INDEX;
     
     int interior_count=0;
     int exterior_count=0;
     int boundary_count=0;
 
     // We start at the same place.
     PHYSBAM_ASSERT( fine_domain.Xmin() == coarse_domain.Xmin() );
     int coarsen_ratio = coarse_domain.Maximum_Edge_Length() / fine_domain.Maximum_Edge_Length();
     
     coarse_voxmap.Fill( EXTERIOR_CELL_TYPE );
     for( RANGE_ITERATOR<d> coarse_iterator( coarse_unpadded_domain ); coarse_iterator.Valid(); coarse_iterator.Next() ){
         int is_full = 0;
         int is_empty = 0;
         
         const T_INDEX& index = coarse_iterator.Index();
         T_INDEX fine_start = (index-1)*coarsen_ratio+1;
         T_INDEX fine_end = (index-1)*coarsen_ratio+T_INDEX::All_Ones_Vector()*coarsen_ratio;
 
         for( RANGE_ITERATOR<d> fine_iterator( RANGE<T_INDEX>(fine_start, fine_end)); fine_iterator.Valid(); fine_iterator.Next()){
                 if( fine_voxmap( fine_iterator.Index() ) == true)
                     is_full++;
                 else
                     is_empty++;                    
             }
             
         coarse_voxmap( index ) = (float)(is_full)/((float)(is_full + is_empty));
     }
 }
 
 
 //#####################################################################
 // Function Multilinear_Interpolation_Stencil
 //#####################################################################
 template<class T,int d>
 STENCIL<T,d> EMBEDDINGTOOLS<T,d>::
 Multilinear_Interpolation_Stencil(const VECTOR<int, d>& cell_index,
                                   const VECTOR<T,d>& multilinear_coordinates)
 {
     typedef STENCIL<T,d> T_STENCIL;
     typedef VECTOR<T,d> TV;
     typedef VECTOR<int,d> T_INDEX;
 
     T_STENCIL interpolation_stencil;
     for(RANGE_ITERATOR<d> iterator(RANGE<T_INDEX>(cell_index,cell_index+1));iterator.Valid();iterator.Next()){
         T_INDEX node_index=iterator.Index();
         T weight=(T)1.;
         for(int v=1;v<=d;v++) weight*=node_index(v)==cell_index(v)?(T)1.-multilinear_coordinates(v):multilinear_coordinates(v);
         interpolation_stencil.Insert(node_index,weight);}
     return interpolation_stencil;
 }
 
 
 //#####################################################################
 // Function Deformation
 //#####################################################################
 template<class T,int d>
 VECTOR<T,d> EMBEDDINGTOOLS<T,d>::
 Deformation(const GRID<VECTOR<T,d> >& domain,
             const VECTOR< ARRAY<T,VECTOR<int,d> >, d>& displacement,
             const VECTOR<int,d>& cell_index,
             const VECTOR<T,d>& multilinear_coordinates,
             bool displace_only)
 {
     typedef STENCIL<T,d> T_STENCIL;
     typedef VECTOR<T,d> TV;
     typedef VECTOR<int,d> T_INDEX;
 
     T_STENCIL interpolation_stencil=Multilinear_Interpolation_Stencil(cell_index,multilinear_coordinates);
     TV result;
     for(int v=1;v<=d;v++) result(v)=interpolation_stencil*displacement(v);
     if( displace_only )
         return result;
     else
         return result+domain.Node(cell_index)+domain.dX*multilinear_coordinates;
 }
 
 //#####################################################################
 // Function Multilinear_Interpolation
 //#####################################################################
 template<class T,int d>
 VECTOR<T,d> EMBEDDINGTOOLS<T,d>::
 Multilinear_Interpolation(const GRID<VECTOR<T,d> >& domain,
                           const VECTOR<int,d>& cell_index,
                           const VECTOR<T,d>& multilinear_coordinates)
 {
     typedef STENCIL<T,d> T_STENCIL;
     typedef VECTOR<T,d> TV;
     typedef VECTOR<int,d> T_INDEX;
 
     return domain.Node(cell_index)+domain.dX*multilinear_coordinates;
 }

template struct EMBEDDINGTOOLS<float, 3>;
