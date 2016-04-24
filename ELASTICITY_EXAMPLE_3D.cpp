//#####################################################################
// Copyright 2014-2016, Michael Doescher, Haixiang Liu, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTICITY_EXAMPLE
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
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
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    tetrahedralized_volume.Update_Number_Nodes();
    
    for(int element=1;element<=mesh.elements.Size();element++){
        
        const VECTOR<int,8>& cell=mesh.elements(element); 
        const int LDB=cell(1);
        const int LDF=cell(2); 
        const int LUB=cell(3);
        const int LUF=cell(4);
        const int RDB=cell(5);
        const int RDF=cell(6);
        const int RUB=cell(7);
        const int RUF=cell(8);

 
            tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LDB,RDB,LUB,LDF));
            tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,RDF,RUF,LDF));
            tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUB,RUB,RUF,RDB));
            tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(LUF,RUF,LDF,LUB));
            tetrahedralized_volume.mesh.elements.Append(VECTOR<int,4>(RDB,LDF,RUF,LUB));}

    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection(particles);
    collection.Add_Structure(&tetrahedralized_volume);

    FILE_UTILITIES::Create_Directory(output_directory+"/"+STRING_UTILITIES::Value_To_String(frame));
    collection.Write(stream_type,output_directory,frame,frame,true);
}
//#####################################################################
template class ELASTICITY_EXAMPLE<VECTOR<float,3> >;
