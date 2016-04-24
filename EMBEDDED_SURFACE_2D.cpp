//#####################################################################
// Copyright 2014-2016, Michael Doescher, Haixiang Liu, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_SURFACE
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include "REGULAR_HYPERCUBE_MESH.h"
#include "EMBEDDED_SURFACE.h"
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include "EMBEDDING_TOOLS.h"
using namespace PhysBAM;    
//#####################################################################
// Function Initialize_Examples
//#####################################################################

template<class TV> void EMBEDDED_SURFACE<TV>::
Load_Embedded_Surface()
{
    if (example==1) // Select test cells to be empty
    {
        // identify full cells
        for (int i=1;i<=mesh.element_coordinates.m;i++) {
            const T_INDEX& cell_index=mesh.element_coordinates(i);
            if (!((T)cell_index(1) > 0.3*(T)size(1)) && ((T)cell_index(1) < 0.7 * (T)size(1)) && ((T)cell_index(2) < 0.7*(T)size(2))) 
                cell_is_full.Set(i);
        }
    }

    if (example==2) { // sphere
        // identify full cells
        for (int i=1;i<=mesh.element_coordinates.m;i++) {
            const T_INDEX& cell_index=mesh.element_coordinates(i);
            cell_is_full.Set(i);
        }

        // generate sphere 
        SPHERE<TV> sphere(TV(size)/(T)2, size.Min()/(T)2.1);
        SEGMENTED_CURVE_2D<T> *segmented_curve = TESSELLATION::Tessellate_Boundary<T>(sphere, 6);
        for (int i=1;i<=segmented_curve->particles.X.m;i++) embedded_vertices.Append(segmented_curve->particles.X(i));
        for (int i=1;i<=segmented_curve->mesh.elements.m;i++) embedded_faces.Append(segmented_curve->mesh.elements(i));
    }
 
    if (example==3) { 
        // generate sphere 
        SPHERE<TV> sphere(TV(size)/(T)2, size.Min()/(T)2.1);
        SEGMENTED_CURVE_2D<T> *segmented_curve = TESSELLATION::Tessellate_Boundary<T>(sphere, 6);
        for (int i=1;i<=segmented_curve->particles.X.m;i++) embedded_vertices.Append(segmented_curve->particles.X(i));
        for (int i=1;i<=segmented_curve->mesh.elements.m;i++) embedded_faces.Append(segmented_curve->mesh.elements(i));

        EMBEDDINGTOOLS<T,d>::Rasterize(embedded_vertices, embedded_faces, dx, mesh, X, cell_is_full);  
        EMBEDDINGTOOLS<T,d>::Generate_Embedding_Map(embedded_vertices, dx, mesh, X, embedding_map);   

        // Select nodes to be kinematically constrained
        for (HASHTABLE_ITERATOR<T_INDEX,int> iterator(node_hash);iterator.Valid();iterator.Next()) {
            T_INDEX node_index=iterator.Key();
            if( node_index(1) < 10 || node_index(1) > size(1)-10 )
                node_is_dirichlet.Set(iterator.Data());
        }
    }


}

//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X)const{

}


//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Write_Output_Files(const int frame) const
{

    // Initialize Particles
    GEOMETRY_PARTICLES<TV> particles;
    particles.array_collection->Add_Elements(X.Size());
    particles.X.Prefix(X.Size())=X;

    // Initialize Domain's Triangulated Surface
    TRIANGULATED_AREA<T>& triangulated_area=*TRIANGULATED_AREA<T>::Create(particles);
    ARRAY<VECTOR<T,4> >subdomain_colors;
    for(int element=1;element<=mesh.elements.Size();element++){
        if (cell_is_full.Contains(element)) {
            const VECTOR<int,4>& cell=mesh.elements(element); 
            triangulated_area.mesh.elements.Append(VECTOR<int,3>(cell(1),cell(3),cell(2)));
            triangulated_area.mesh.elements.Append(VECTOR<int,3>(cell(2),cell(3),cell(4)));
            if( subdomain_membership.Get( element ) % 2 == 0){
                subdomain_colors.Append(VECTOR<T,4>(1.0,0.0,0.0,1.0));
                subdomain_colors.Append(VECTOR<T,4>(1.0,0.0,0.0,1.0));
            }
            else{
                subdomain_colors.Append(VECTOR<T,4>(0.0,0.0,1.0,1.0));
                subdomain_colors.Append(VECTOR<T,4>(0.0,0.0,1.0,1.0));            
            }
        }
    }

    // Initialize the segmented curve
    SEGMENTED_CURVE_2D<T>& embedded_segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);
    if (example==2) {
    int particle_offset=particles.X.m;
    for (int i=1;i<=embedded_vertices.m;i++) {
        int p=particles.array_collection->Add_Element();
        particles.X(p)=embedded_vertices(i);
    }
    for (int i=1;i<=embedded_faces.m;i++) embedded_segmented_curve.mesh.elements.Append(embedded_faces(i)+particle_offset);
    }

    if (example==3) {
    // Initialize the segmented curve
        int particle_offset=particles.X.m;
        for (int i=1;i<=embedding_map.m;i++) {
            TV multilinear_coordinates=embedding_map(i).y;
            int cell=embedding_map(i).x;
            TV min_cell_node=X(mesh.elements(cell)(1));
            TV max_cell_node=X(mesh.elements(cell)(4));
            TV world_space_coordinates=min_cell_node+(multilinear_coordinates*(max_cell_node-min_cell_node));
            int p=particles.array_collection->Add_Element();
            particles.X(p)=(world_space_coordinates);
        }
        for (int i=1;i<=embedded_faces.m;i++) embedded_segmented_curve.mesh.elements.Append(embedded_faces(i)+particle_offset);
    }

    // Initialize Nodes
    FREE_PARTICLES<TV>& exterior_nodes=*FREE_PARTICLES<TV>::Create(particles);
    for(HASHTABLE_ITERATOR<int> iterator(node_is_exterior); iterator.Valid(); iterator.Next()){
        exterior_nodes.nodes.Append( iterator.Key() );
    }

    // Initialize Dirichlet nodes
    FREE_PARTICLES<TV>& dirichlet_nodes=*FREE_PARTICLES<TV>::Create(particles);
    for(HASHTABLE_ITERATOR<int> iterator(node_is_dirichlet); iterator.Valid(); iterator.Next()){
        dirichlet_nodes.nodes.Append( iterator.Key() );
    }


    // Update_Number_Nodes
    triangulated_area.Update_Number_Nodes();    
    embedded_segmented_curve.Update_Number_Nodes();
    exterior_nodes.Update_Number_Nodes();
    dirichlet_nodes.Update_Number_Nodes();


    // Add everything to collection for output
    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(particles);
    collection.Add_Structure(&triangulated_area);
    collection.Add_Structure(&exterior_nodes);
    collection.Add_Structure(&dirichlet_nodes);
    collection.Add_Structure(&embedded_segmented_curve);

    FILE_UTILITIES::Create_Directory(output_directory+"/"+STRING_UTILITIES::Value_To_String(frame));
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+STRING_UTILITIES::Value_To_String(frame)+"/stress_map_of_triangulated_area_1", subdomain_colors);
    collection.Write(stream_type,output_directory,frame,frame,true);
}

//#####################################################################
template class EMBEDDED_SURFACE<VECTOR<float,2> >;

