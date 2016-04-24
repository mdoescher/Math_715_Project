//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COROTATED
//#####################################################################
#ifndef __COROTATED__
#define __COROTATED__

#include "MATERIAL_MODEL.h"

namespace PhysBAM{

template<class T,int d>
class COROTATED
{
protected:

    typedef T T_SCALAR;
    enum {dim=d};
    enum {is_mooney_rivlin=0};

public:

    // Consider making protected
    static T M(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return (Sigma-1.).Trace();}

    static DIAGONAL_MATRIX<T,d> Q_hat(const DIAGONAL_MATRIX<T,d>& Sigma)
    {return DIAGONAL_MATRIX<T,d>::Identity_Matrix();}

    static DIAGONAL_MATRIX<T,d> Deviatoric_P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu)
    {return (Sigma-1.)*(2.*mu);}

    static ROTATED_STRESS_DERIVATIVE<T,2> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T mu,const T kappa)
    {
        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;

        rotated_dPdF.a1111=2.*mu+kappa;
        rotated_dPdF.a2222=2.*mu+kappa;
        rotated_dPdF.a1122=kappa;

        rotated_dPdF.a1212=2.*mu+(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(1,1)+Sigma(2,2));
        rotated_dPdF.a1221=-(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(1,1)+Sigma(2,2));

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T mu,const T kappa)
    {
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;

        rotated_dPdF.a1111=2.*mu+kappa;
        rotated_dPdF.a2222=2.*mu+kappa;
        rotated_dPdF.a3333=2.*mu+kappa;
        rotated_dPdF.a1122=kappa;
        rotated_dPdF.a1133=kappa;
        rotated_dPdF.a2233=kappa;

        rotated_dPdF.a1212=2.*mu+(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(1,1)+Sigma(2,2));
        rotated_dPdF.a1221=-(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(1,1)+Sigma(2,2));
        rotated_dPdF.a1313=2.*mu+(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(1,1)+Sigma(3,3));
        rotated_dPdF.a1331=-(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(1,1)+Sigma(3,3));
        rotated_dPdF.a2323=2.*mu+(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(2,2)+Sigma(3,3));
        rotated_dPdF.a2332=-(kappa*(Sigma-1.).Trace()-2.*mu)/(Sigma(2,2)+Sigma(3,3));

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,2> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,2>& Sigma,const T p,const T mu,const T alpha)
    {
        ROTATED_STRESS_DERIVATIVE<T,2> rotated_dPdF;

        rotated_dPdF.a1111=2.*mu;
        rotated_dPdF.a2222=2.*mu;
        rotated_dPdF.a1122=T();

        rotated_dPdF.a1212=2.*mu+(alpha*p-2.*mu)/(Sigma(1,1)+Sigma(2,2));
        rotated_dPdF.a1221=-(alpha*p-2.*mu)/(Sigma(1,1)+Sigma(2,2));

        return rotated_dPdF;
    }

    static ROTATED_STRESS_DERIVATIVE<T,3> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T p,const T mu,const T alpha)
    {
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;

        rotated_dPdF.a1111=2.*mu;
        rotated_dPdF.a2222=2.*mu;
        rotated_dPdF.a3333=2.*mu;
        rotated_dPdF.a1122=T();
        rotated_dPdF.a1133=T();
        rotated_dPdF.a2233=T();

        rotated_dPdF.a1212=2.*mu+(alpha*p-2.*mu)/(Sigma(1,1)+Sigma(2,2));
        rotated_dPdF.a1221=-(alpha*p-2.*mu)/(Sigma(1,1)+Sigma(2,2));
        rotated_dPdF.a1313=2.*mu+(alpha*p-2.*mu)/(Sigma(1,1)+Sigma(3,3));
        rotated_dPdF.a1331=-(alpha*p-2.*mu)/(Sigma(1,1)+Sigma(3,3));
        rotated_dPdF.a2323=2.*mu+(alpha*p-2.*mu)/(Sigma(2,2)+Sigma(3,3));
        rotated_dPdF.a2332=-(alpha*p-2.*mu)/(Sigma(2,2)+Sigma(3,3));

        return rotated_dPdF;
    }

//#################################################
    static DIAGONAL_MATRIX<T,d> Deviatoric_P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu_10,const T mu_01)
    {PHYSBAM_FATAL_ERROR();}
    static ROTATED_STRESS_DERIVATIVE<T,d> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T mu_10,const T mu_01,const T kappa) {PHYSBAM_FATAL_ERROR();}
    static ROTATED_STRESS_DERIVATIVE<T,d> Augmented_Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu_10,const T mu_01,const T alpha) {PHYSBAM_FATAL_ERROR();}
//#################################################

};

}
#endif
