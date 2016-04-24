#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/BOX_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>

#include "EMBEDDING_TOOLS.h"


using namespace PhysBAM;


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

        // traverse the list of segments and identify which cells are intersected by each segment
        for(int t=1;t<=faces.m;t++){
            RANGE<TV> simplex_box=RANGE<TV>::Empty_Box();
            for (int v=1;v<=d;v++) simplex_box.Enlarge_To_Include_Point(vertices(faces(t)(v)));
            ARRAY<int> intersection_list;
            box_hierarchy.Intersection_List(simplex_box, intersection_list);
            SEGMENT_2D<T> segment(VECTOR<T,2>(vertices(faces(t)(1)).x,vertices(faces(t)(1)).y ),VECTOR<T,2>(vertices(faces(t)(2)).x,vertices(faces(t)(2)).y ));
            for (int i=1;i<=intersection_list.m;i++) {
                RANGE<TV> cell = boxes(intersection_list(i));   
                if(INTERSECTION::Intersects<T>(cell,segment,1e-5)) cell_is_full.Set(intersection_list(i));
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


template struct EMBEDDINGTOOLS<float, 2>;
