//#####################################################################
// Copyright 2014-2016, Nathan Mitchell, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CG_VECTOR
//#####################################################################
#ifndef __CG_VECTOR__
#define __CG_VECTOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>

using namespace PhysBAM;

namespace PhysBAM{

template<class TV>
class CG_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::ELEMENT>
{
    typedef typename TV::ELEMENT T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;

    ARRAY<TV>& array;
public:
    CG_VECTOR(ARRAY<TV>& array_input) : array(array_input) {}

    static ARRAY<TV>& Array(BASE& base_array)
    {return ((CG_VECTOR&)(base_array)).array;}

    static const ARRAY<TV>& Array(const BASE& base_array)
    {return ((const CG_VECTOR&)(base_array)).array;}

    BASE& operator+=(const BASE& bv)
    {array+=Array(bv);return *this;}

    BASE& operator-=(const BASE& bv)
    {array-=Array(bv);return *this;}

    BASE& operator*=(const T a)
    {array*=a;return *this;}

    void Copy(const T c,const BASE& bv)
    {ARRAY<TV>::Copy(c,Array(bv),array);}

    void Copy(const T c1,const BASE& bv1,const BASE& bv2)
    {ARRAY<TV>::Copy(c1,Array(bv1),Array(bv2),array);}

    int Raw_Size() const
    {return array.Flattened().m;}
    
    T& Raw_Get(int i)
    {return array.Flattened()(i);}

//#####################################################################
//#####################################################################
};
}
#endif
