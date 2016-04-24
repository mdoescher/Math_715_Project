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
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include "REGULAR_HYPERCUBE_MESH.h"
#include "EMBEDDED_SURFACE.h"
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include "EMBEDDING_TOOLS.h"

using namespace PhysBAM;   

//#####################################################################
// Function Rescale_Model_And_Domain
//#####################################################################
template<class TV> void Rescale_Model_And_Domain(ARRAY<TV>& embedded_vertices, RANGE<TV>& domain, RANGE<TV>& model_box, VECTOR<int,3>& size, const typename TV::SCALAR dx)
{
    // model bounding_box
    model_box=RANGE<TV>::Empty_Box();
    for (int i=1;i<=embedded_vertices.m;i++) model_box.Enlarge_To_Include_Point(embedded_vertices(i));
    TV model_size=model_box.Edge_Lengths(); 

    // resize and center domain
    domain=RANGE<TV>(TV(),TV::All_Ones_Vector()/(model_size.Min()/model_size));
    domain=RANGE<TV>(TV()-domain.Center(), domain.Maximum_Corner()-domain.Center());

    // resize model
    for (int i=1;i<=embedded_vertices.m;i++) {
        embedded_vertices(i)-=model_box.Center();
        embedded_vertices(i)/=model_size.Min();
        embedded_vertices(i)*=0.999;}

    // resize the bounding box of the model
    TV corner=model_box.Minimum_Corner();
    corner-=model_box.Center();
    corner/=model_size.Min();
    corner*=0.999;
    model_box=RANGE<TV>(corner,-corner);

    // update size
    TV t_size=ceil(domain.Edge_Lengths()/dx);
    size=VECTOR<int,3>((int)t_size.x, (int)t_size.y, (int)t_size.z);

}
//#####################################################################
// Function Load_Embedded_Surface
// #####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Load_Embedded_Surface()
{
    if (example==1 || example == 2) { // sphere uisng PhysBAM SPHERE
        SPHERE<TV> sphere(TV(), (T)1.);
        TRIANGULATED_SURFACE<T> *embedded_triangulated_surface = TESSELLATION::Tessellate_Boundary<T>(sphere, 4);
        for (int i=1;i<=embedded_triangulated_surface->particles.X.m;i++) embedded_vertices.Append(embedded_triangulated_surface->particles.X(i));
        for (int i=1;i<=embedded_triangulated_surface->mesh.elements.m;i++) embedded_faces.Append(embedded_triangulated_surface->mesh.elements(i));
    }

    if (example==3 || example == 4) { // load sphere model from file 
        std::string input_filename="/store1/sifakis_group/sphere.tri";
        TRIANGULATED_SURFACE<T>& embedded_triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,input_filename,embedded_triangulated_surface);
        for (int i=1;i<=embedded_triangulated_surface.particles.X.m;i++) embedded_vertices.Append(embedded_triangulated_surface.particles.X(i));
        for (int i=1;i<=embedded_triangulated_surface.mesh.elements.m;i++) embedded_faces.Append(embedded_triangulated_surface.mesh.elements(i));
    }

    if (example==5) { // lion2
        std::string input_filename="models/lion2.tri";
        TRIANGULATED_SURFACE<T>& embedded_triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,input_filename,embedded_triangulated_surface);
        for (int i=1;i<=embedded_triangulated_surface.particles.X.m;i++) embedded_vertices.Append(embedded_triangulated_surface.particles.X(i));
        for (int i=1;i<=embedded_triangulated_surface.mesh.elements.m;i++) embedded_faces.Append(embedded_triangulated_surface.mesh.elements(i));
        ROTATION<TV> rotation(T(3.14159/2.),TV(1,0,0));
        for (int i=1;i<=embedded_vertices.m;i++) embedded_vertices(i)=rotation.Rotate(embedded_vertices(i));
    }

    if (example==6) { //xyzrgb_dragon
        std::string input_filename="models/xyzrgb_dragon_7mil.tri";
        TRIANGULATED_SURFACE<T>& embedded_triangulated_surface=*TRIANGULATED_SURFACE<T>::Create();
        FILE_UTILITIES::Read_From_File(stream_type,input_filename,embedded_triangulated_surface);
        for (int i=1;i<=embedded_triangulated_surface.particles.X.m;i++) embedded_vertices.Append(embedded_triangulated_surface.particles.X(i));
        for (int i=1;i<=embedded_triangulated_surface.mesh.elements.m;i++) embedded_faces.Append(embedded_triangulated_surface.mesh.elements(i));
        ROTATION<TV> rotation(T(3.14159/2.),TV(1,0,0));
        ROTATION<TV> rotation2(T(3.14159),TV(0,1,0));
        for (int i=1;i<=embedded_vertices.m;i++) {
            embedded_vertices(i)=rotation.Rotate(embedded_vertices(i));
            embedded_vertices(i)=rotation2.Rotate(embedded_vertices(i));}
    }

    if (example==11) {  // two spheres one mildly perturbed the other wildly perturbed
        SPHERE<TV> sphere1(TV(-1.,0.,0.), .8);
        SPHERE<TV> sphere2(TV( 1.,0.,0.), .8);
        TRIANGULATED_SURFACE<T> *embedded_triangulated_surface1 = TESSELLATION::Tessellate_Boundary<T>(sphere1,6);
        TRIANGULATED_SURFACE<T> *embedded_triangulated_surface2 = TESSELLATION::Tessellate_Boundary<T>(sphere2,6);
        for (int i=1;i<=embedded_triangulated_surface1->particles.X.m;i++)   embedded_vertices.Append(embedded_triangulated_surface1->particles.X(i));
        for (int i=1;i<=embedded_triangulated_surface1->mesh.elements.m;i++) embedded_faces.Append(embedded_triangulated_surface1->mesh.elements(i));
        int offset=embedded_vertices.m;
        for (int i=1;i<=embedded_triangulated_surface2->particles.X.m;i++)   embedded_vertices.Append(embedded_triangulated_surface2->particles.X(i));
        for (int i=1;i<=embedded_triangulated_surface2->mesh.elements.m;i++) {
            VECTOR<int,3> face=embedded_triangulated_surface2->mesh.elements(i);
            embedded_faces.Append(face+VECTOR<int,3>(offset,offset,offset));}
    }



    int i=model_ranges.Append(RANGE<TV>());
    Rescale_Model_And_Domain(embedded_vertices, domain, model_ranges(i), size, dx);
}
//#####################################################################
// Function Select_Kinematic_Nodes
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Select_Kinematic_Nodes()
{   // this uses an array of axis aligned ranges to select regions for scripting : model_ranges(1) is the entire model 
    kinematic_nodes.Resize(2);
    if (example >= 1 && example <=4) { // the spheres
        for (int i=1;i<=X.m;i++) {
            if (X(i)(1)<=model_ranges(1).Minimum_Corner()(1)) {
                node_is_dirichlet.Set(i);
                kinematic_nodes(1).Append(i); // left side
            }
            if (X(i)(1)>model_ranges(1).Maximum_Corner()(1)) {
                node_is_dirichlet.Set(i);
                kinematic_nodes(2).Append(i);
            }
        }
    }

    if (example==5){
        kinematic_nodes.Resize(2);
        // feet
        for (int i=1;i<=X.m;i++) {
            if ((X(i)(2)<=model_ranges(1).Minimum_Corner()(2))) {
                node_is_dirichlet.Set(i);
                kinematic_nodes(1).Append(i);}}

        // head
        T left   = 0.0;      T right = 0.24;
        T bottom = 0.55;     T top   = 1.;
        T back   = 0.3;      T front = 0.8;
        const TV sm_el = model_ranges(1).Edge_Lengths();
        TV head_min_corner=model_ranges(1).Minimum_Corner() + TV(sm_el.x*left,  sm_el.y*bottom, sm_el.z*back);
        TV head_max_corner=model_ranges(1).Minimum_Corner() + TV(sm_el.x*right, sm_el.y*top,    sm_el.z*front);
        RANGE<TV> head(head_min_corner, head_max_corner);
        int j=model_ranges.Append(head);
        for (int i=1;i<=X.m;i++) {
            if (model_ranges(j).Lazy_Inside(X(i))) {
                node_is_dirichlet.Set(i);
                kinematic_nodes(2).Append(i);
            }
        }
    }

    if (example==6) {
        kinematic_nodes.Resize(2);

        // head
        T left   = 0.05;      T right = 0.25;
        T bottom = 0.55;     T top   = .85;
        T back   = 0.35;      T front = 0.55;
        const TV sm_el = model_ranges(1).Edge_Lengths();
        TV head_min_corner=model_ranges(1).Minimum_Corner() + TV(sm_el.x*left,  sm_el.y*bottom, sm_el.z*back);
        TV head_max_corner=model_ranges(1).Minimum_Corner() + TV(sm_el.x*right, sm_el.y*top,    sm_el.z*front);
        RANGE<TV> head(head_min_corner, head_max_corner);
        int j=model_ranges.Append(head);
        for (int i=1;i<=X.m;i++) {
            if (model_ranges(j).Lazy_Inside(X(i))) {
                node_is_dirichlet.Set(i);
                kinematic_nodes(2).Append(i);
            }
        }
    }

    if (example==11) {  // two spheres one mildly perturbed the other wildly perturbed
        kinematic_nodes.Resize(1);
        SPHERE<TV> sphere1(TV(-.5,0.,0.), .1); // just anchor the center of the spheres in position
        SPHERE<TV> sphere2(TV(.5,0.,0.), .1); // just anchor the spheres in position
        for (int i=1;i<=X.m;i++) {
                if (sphere1.Lazy_Inside(X(i))) {node_is_dirichlet.Set(i);}
                if (sphere2.Lazy_Inside(X(i))) {
                    node_is_dirichlet.Set(i);
                    kinematic_nodes(1).Append(i);
                }
        }
        RANDOM_NUMBERS<T> random(1);
        SPHERE<TV> sphere3(TV(-.5,0.,0.),sphere1.radius+.001);
        for (int i=1;i<=X.m;i++) {
            if (!(node_is_dirichlet.Contains(i))){
                if (X(i)(1)<0.) {
                    if (!sphere3.Lazy_Inside(X(i)) && random.Get_Uniform_Number(0.,1.)<0.05){
                    X(i)+=(X(i)-sphere1.center)*perturbation_magnitude; }
                }
            }
        }
    }




}

//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class TV> void EMBEDDED_SURFACE<TV>::
Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X)const{
    // stretch 
    if (example== 1) { 
        for (int j=1;j<=kinematic_nodes(1).m;j++) {
            const int  node = kinematic_nodes(1)(j);
            X(node)(1)-=time*dx*T(0.3);
        }
        for (int j=1;j<=kinematic_nodes(2).m;j++) {
            const int  node = kinematic_nodes(2)(j);
            X(node)(1)+=time*dx*T(0.3);
        }
    }

    // rotate ends
    if (example==2 || example==4) {
        T theta=(T)3.14159/(T)40./(T)20.;
        MATRIX<T,3> M=MATRIX<T,3>(1, 0, 0, 0, cos(theta), T(-1)*sin(theta), 0, sin(theta), cos(theta));
        for (int j=1;j<=kinematic_nodes(1).m;j++) {
            const int  node = kinematic_nodes(1)(j);
            X(node)=M*X(node);
        }
        theta*=(T)-1.;
        M=MATRIX<T,3>(1, 0, 0, 0, cos(theta), T(-1)*sin(theta), 0, sin(theta), cos(theta));
        for (int j=1;j<=kinematic_nodes(2).m;j++) {
            const int  node = kinematic_nodes(2)(j);
            X(node)=M*X(node);
        }
    }
    if (example==3) { 
        for (int j=1;j<=kinematic_nodes(1).m;j++) {
            const int  node = kinematic_nodes(1)(j);
            X(node)(1)-=time*dx*T(0.05);
        }
        for (int j=1;j<=kinematic_nodes(2).m;j++) {
            const int  node = kinematic_nodes(2)(j);
            X(node)(1)+=time*dx*T(0.05);
        }
    }
  

    // lion 
    // model_ranges(1)=whole model, (2)=head
    // kinematic_nodes(1)=ground plane (2)=head
    // macro_block_elasticity_nocona -example 5 -dx .1 -constant_mu 10 -constant_lambda 90 -o embedded_surface_lion -last_frame 360 -newton_iterations 1 -newton_tolerance 1e-4 -CG_iterations 20 
    if (example==5) {
        T sign=1.;if(time > 5.) sign=-1.;
        T max_rotation_angle=3.14159*5./12.;
        T theta = max_rotation_angle/5./frame_rate*sign;
        ROTATION<TV> rotation(theta,TV::Axis_Vector(1));
        for (int i=1;i<=kinematic_nodes(2).m;i++) {
           int node=kinematic_nodes(2)(i);
            X(node)-=model_ranges(2).Center();
            X(node)=rotation.Rotate(X(node));
            X(node)+=model_ranges(2).Center();
        }
    }

    // dragon
    // model_ranges(1)=whole model, (2)=head
    // kinematic_nodes(1) empty, (2)=head
    // macro_block_elasticity_nocona -example 6 -dx .05 -constant_mu 10 -constant_lambda 90 -o embedded_surface_dragon -last_frame 120 -newton_iterations 1 -newton_tolerance 1e-4 -CG_iterations 20 
    if (example==6) {
        //Parabolic Helix
        T y_radius=.03;
        T z_radius=.04;
        T forward_rate=-.002;
        TV h(forward_rate*time,y_radius*cos(time),z_radius*sin(time));
        for (int i=1;i<=kinematic_nodes(2).m;i++) {
           int node=kinematic_nodes(2)(i);
            X(node)+=h;
        }
    }

    if (example==11) {    // two spheres 
        for (int i=1;i<=kinematic_nodes(1).m;i++) {
            X(kinematic_nodes(1)(i))+=TV::Axis_Vector(1)*rate;
        }
    }
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

    // Tetrahedralized Volume
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    tetrahedralized_volume.Update_Number_Nodes();

    for(int element=1;element<=mesh.elements.Size();element++){
        if (cell_is_full.Contains(element)) {
            const VECTOR<int,8>& cell=mesh.elements(element); 
            const int LDB=cell(1);
            const int LDF=cell(2); 
            const int LUB=cell(3);
            const int LUF=cell(4);
            const int RDB=cell(5);
            const int RDF=cell(6);
            const int RUB=cell(7);
            const int RUF=cell(8);

            const VECTOR<int,3>& node=node_to_ijk.Get(cell(1));
            int i,j,ij;node.Get(i,j,ij);

            if((i+j+ij)%2==0){
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,LUB,LDF));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,RDF,RUF,LDF));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUB,RUB,RUF,RDB));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,LDF,LUB));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,LDF,RUF,LUB));}
            else{
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,RUB,RDF));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,LUB,LUF,RUB));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RDF,LDF,LDB));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,RDF,RUB));
                tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,LDB,RUB,RDF));}}
    }

    // Embedded Triangulated Surface 
    TRIANGULATED_SURFACE<T>& embedded_triangulated_surface=*TRIANGULATED_SURFACE<T>::Create(particles);
    int particle_offset=particles.X.m;
    ARRAY<VECTOR<int,d> > offset_faces(embedded_faces.m);
    for (int i=1;i<=embedded_faces.m;i++) offset_faces(i)=embedded_faces(i)+particle_offset;
    embedded_triangulated_surface.mesh.Initialize_Mesh(offset_faces.m,offset_faces);
    for (int i=1;i<=embedding_map.m;i++) {
        TV coordinates=embedding_map(i).y;
        TV deformed_position;
        int node=1;
        for( RANGE_ITERATOR<d> node_iterator( RANGE<T_INDEX>::Unit_Box() ); node_iterator.Valid(); node_iterator.Next(), node++ ){
            T per_node_weight=1;
            for( int w=1;w<=d;w++ )
                if( node_iterator.Index()(w)) per_node_weight *= coordinates(w);
                else per_node_weight *= 1.0-coordinates(w);
            deformed_position += X(mesh.elements(embedding_map(i).x)(node)) * per_node_weight;}
        int p=particles.array_collection->Add_Element();
        particles.X(p)=deformed_position;
    } 


    // dirichlet Points
    FREE_PARTICLES<TV>& dirichlet_points=*FREE_PARTICLES<TV>::Create(particles);
    for (HASHTABLE_ITERATOR<int> iterator(node_is_dirichlet);iterator.Valid();iterator.Next()) {
        const int node=iterator.Key();
        const int p=particles.array_collection->Add_Element();
        particles.X(p)=X(node);
        dirichlet_points.nodes.Append(p);
    }

    // update number of nodes
    tetrahedralized_volume.Update_Number_Nodes();
    embedded_triangulated_surface.Update_Number_Nodes();
    dirichlet_points.Update_Number_Nodes();

    // Add everything to collection for output
    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(particles);
    if (tetrahedralized_volume.mesh.elements.Size()!=0) collection.Add_Structure(&tetrahedralized_volume);
    if (embedded_triangulated_surface.mesh.elements.Size()!=0) collection.Add_Structure(&embedded_triangulated_surface);
    if (dirichlet_points.nodes.m!=0) collection.Add_Structure(&dirichlet_points);

    FILE_UTILITIES::Create_Directory(output_directory+"/"+STRING_UTILITIES::Value_To_String(frame));
    collection.Write(stream_type,output_directory,frame,frame,true);
}

//#####################################################################
template class EMBEDDED_SURFACE<VECTOR<float,3> >;

