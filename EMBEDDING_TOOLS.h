//#####################################################################
// Copyright 2010-2013, Michael Doescher, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EMBEDDING_TOOLS_H__
#define __EMBEDDING_TOOLS_H__

#include <iostream>
#include <iomanip>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include "CELL_TYPE.h"
#include "RANGE_ITERATOR.h"
#include "STENCIL.h"
#include "REGULAR_HYPERCUBE_MESH.h"

namespace PhysBAM {

    template< class T, int d>
    struct EMBEDDINGTOOLS{

//#####################################################################
// Function Get_Grid_Bounds
//#####################################################################
        
        static void Get_Grid_Bounds( const ARRAY<VECTOR<T,d> >& vertices, const ARRAY<VECTOR<int,d> >& triangles, const T dx, VECTOR<int, d>& cell_bounds, VECTOR<T, d>& min_corner);
        
        
//#####################################################################
// Function Rasterize
//#####################################################################
        
        static void Rasterize( const ARRAY<VECTOR<T,d> >& vertices, const ARRAY<VECTOR<int,d> >& triangles,
                               const GRID<VECTOR<T,d> >& domain, const RANGE<VECTOR<int,d> >& unpadded_domain,
                               const RANGE<VECTOR<int,d> >& padded_domain, ARRAY< bool, VECTOR<int, d> >& voxmap,
                               ARRAY< bool, VECTOR<int, d> >& voxmap_node);
        
        
        static void Rasterize(const ARRAY<VECTOR<T,d> >& embedded_vertices, const ARRAY<VECTOR<int,d> >& embedded_faces, const T dx, 
                              const REGULAR_HYPERCUBE_MESH<d>& mesh, const ARRAY<VECTOR<T,d> >& X, HASHTABLE<int>& cell_is_full);


//#####################################################################
// Function Generate_Embedding_Map
//#####################################################################
        
        void static Generate_Embedding_Map(const ARRAY<VECTOR<T,d> > &vertices, const T dx, const REGULAR_HYPERCUBE_MESH<d>& mesh, 
                               const ARRAY<VECTOR<T,d> >& X, ARRAY<PAIR<int,VECTOR<T,d> > >& embedding_map);


//#####################################################################
// Function Coarsen
//#####################################################################
        
        static void Coarsen( const GRID<VECTOR<T,d> >& fine_domain, const RANGE<VECTOR<int,d> >& fine_unpadded_domain,
                             const RANGE<VECTOR<int,d> >& fine_padded_domain, const ARRAY< bool, VECTOR<int, d> >& fine_voxmap,
                             const GRID<VECTOR<T,d> >& coarse_domain, const RANGE<VECTOR<int,d> >& coarse_unpadded_domain,
                             const RANGE<VECTOR<int,d> >& coarse_padded_domain, ARRAY< CELL_TYPE , VECTOR<int, d> >& coarse_voxmap);
        
//#####################################################################
// Function Coarsen_Density
//#####################################################################
        
        static void Coarsen_Density( const GRID<VECTOR<T,d> >& fine_domain, const RANGE<VECTOR<int,d> >& fine_unpadded_domain,
                                     const RANGE<VECTOR<int,d> >& fine_padded_domain, const ARRAY< bool, VECTOR<int, d> >& fine_voxmap,
                                     const GRID<VECTOR<T,d> >& coarse_domain, const RANGE<VECTOR<int,d> >& coarse_unpadded_domain,
                                     const RANGE<VECTOR<int,d> >& coarse_padded_domain, ARRAY< float , VECTOR<int, d> >& coarse_voxmap);
        
//#####################################################################
// Function Deformation
//#####################################################################
        
        static VECTOR<T,d> Deformation(const GRID<VECTOR<T,d> >& domain,
                                       const VECTOR< ARRAY< T, VECTOR<int,d> >, d >& displacement,
                                       const VECTOR<int,d>& cell_index,
                                       const VECTOR<T,d>& multilinear_coordinates,
                                       bool displace_only);

//#####################################################################
// Function Multilinear_Interpolation
//#####################################################################
        
        static VECTOR<T,d> Multilinear_Interpolation(const GRID<VECTOR<T,d> >& domain,
                                                     const VECTOR<int,d>& cell_index,
                                                     const VECTOR<T,d>& multilinear_coordinates);
        
//#####################################################################
// Function Multilinear_Interpolation_Stencil
//#####################################################################
        static STENCIL<T,d> Multilinear_Interpolation_Stencil(const VECTOR<int,d>& cell_index,
                                                              const VECTOR<T,d>& multilinear_coordinates);
        
//#####################################################################
    };
}
#endif
