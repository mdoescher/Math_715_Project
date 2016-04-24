//#####################################################################
// Copyright 2014-2016, Haixiang Liu, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTICITY_EXAMPLE
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "REGULAR_HYPERCUBE_MESH.h"
#include "ELASTICITY_EXAMPLE.h"
using namespace PhysBAM; 
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void ELASTICITY_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    GEOMETRY_PARTICLES<TV> particles;
    particles.array_collection->Add_Elements(X.Size());
    particles.X.Prefix(X.Size())=X;
    TRIANGULATED_AREA<T>& triangulated_area=*TRIANGULATED_AREA<T>::Create(particles);
    for(int element=1;element<=mesh.elements.Size();element++){
        const VECTOR<int,4>& cell=mesh.elements(element); 
        triangulated_area.mesh.elements.Append(VECTOR<int,3>(cell(1),cell(3),cell(2)));
        triangulated_area.mesh.elements.Append(VECTOR<int,3>(cell(2),cell(3),cell(4)));}

    triangulated_area.Update_Number_Nodes();

    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(particles);
    collection.Add_Structure(&triangulated_area);

    FILE_UTILITIES::Create_Directory(output_directory+"/"+STRING_UTILITIES::Value_To_String(frame));
    collection.Write(stream_type,output_directory,frame,frame,true);
}
//#####################################################################
template class ELASTICITY_EXAMPLE<VECTOR<float,2> >;
