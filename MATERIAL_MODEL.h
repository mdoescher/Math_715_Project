//#####################################################################
// Copyright 2011, Taylor Patterson, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATERIAL_MODEL
//#####################################################################
#ifndef __MATERIAL_MODEL__
#define __MATERIAL_MODEL__

#include "ROTATED_STRESS_DERIVATIVE.h"

namespace PhysBAM{

template<class T_MATERIAL>
class MATERIAL_MODEL:public T_MATERIAL
{
    typedef typename T_MATERIAL::T_SCALAR T;
    enum {d=T_MATERIAL::dim};

public:

    static bool Is_Mooney_Rivlin()
    {return T_MATERIAL::is_mooney_rivlin;}

    // Non Mooney-Rivlin
    static DIAGONAL_MATRIX<T,d> P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu,const T alpha)
    {return Deviatoric_P_hat(Sigma,mu)+alpha*p*Q_hat(Sigma);}

    // Mooney-Rivlin
    static DIAGONAL_MATRIX<T,d> P_hat(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu_10,const T mu_01,const T alpha)
    {return Deviatoric_P_hat(Sigma,mu_10,mu_01)+alpha*p*Q_hat(Sigma);}

    static T q(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T alpha,const T alpha_squared_over_kappa)
    {return -alpha*M(Sigma)+alpha_squared_over_kappa*p;}

    static MATRIX<T,d> dP_hat(const ROTATED_STRESS_DERIVATIVE<T,d>& rotated_dPdF,const DIAGONAL_MATRIX<T,d>& Q_hat,const MATRIX<T,d>& dF_hat,const T dp,const T alpha)
    {return rotated_dPdF.dP_hat(dF_hat)+alpha*dp*Q_hat;}

    static T dq(const DIAGONAL_MATRIX<T,d>& Q_hat,const MATRIX<T,d>& dF_hat,const T dp,const T alpha,const T alpha_squared_over_kappa)
    {return -alpha*Q_hat.Times_Transpose(dF_hat).Trace()+alpha_squared_over_kappa*dp;}

    // Non Mooney-Rivlin
    static ROTATED_STRESS_DERIVATIVE<T,d> Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu,const T kappa,const T alpha,const bool apply_definiteness_fix)
    {
        PHYSBAM_ASSERT(!Is_Mooney_Rivlin());

        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF=Augmented_Rotated_Stress_Derivative(Sigma,p,mu,alpha);

#if 1
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_indefinite=Rotated_Stress_Derivative(Sigma,mu,std::min((T)5.*mu,kappa));
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_incremental_fix=rotated_dPdF_indefinite;
        rotated_dPdF_incremental_fix.Make_Positive_Definite();
        rotated_dPdF_incremental_fix-=rotated_dPdF_indefinite;
        if(apply_definiteness_fix) rotated_dPdF+=rotated_dPdF_incremental_fix;
#elif 1
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_indefinite=Rotated_Stress_Derivative(Sigma,mu,kappa);
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_incremental_fix=rotated_dPdF_indefinite;
        rotated_dPdF_incremental_fix.Make_Positive_Definite();
        rotated_dPdF_incremental_fix-=rotated_dPdF_indefinite;
        if(apply_definiteness_fix) rotated_dPdF+=rotated_dPdF_incremental_fix;
#elif 0
        if(apply_definiteness_fix) rotated_dPdF.Make_Positive_Definite();
#endif

        return rotated_dPdF;
    }

    // Mooney-Rivlin
    static ROTATED_STRESS_DERIVATIVE<T,d> Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& Sigma,const T p,const T mu_10,const T mu_01,const T kappa,const T alpha,const bool apply_definiteness_fix)
    {
        PHYSBAM_ASSERT(Is_Mooney_Rivlin());

        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF=Augmented_Rotated_Stress_Derivative(Sigma,p,mu_10,mu_01,alpha);

#if 1
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_indefinite=Rotated_Stress_Derivative(Sigma,mu_10,mu_01,std::min((T)10.*(mu_10+mu_01),kappa));
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_incremental_fix=rotated_dPdF_indefinite;
        rotated_dPdF_incremental_fix.Make_Positive_Definite();
        rotated_dPdF_incremental_fix-=rotated_dPdF_indefinite;
        if(apply_definiteness_fix) rotated_dPdF+=rotated_dPdF_incremental_fix;
#elif 1
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_indefinite=Rotated_Stress_Derivative(Sigma,mu_10,mu_01,kappa);
        ROTATED_STRESS_DERIVATIVE<T,d> rotated_dPdF_incremental_fix=rotated_dPdF_indefinite;
        rotated_dPdF_incremental_fix.Make_Positive_Definite();
        rotated_dPdF_incremental_fix-=rotated_dPdF_indefinite;
        if(apply_definiteness_fix) rotated_dPdF+=rotated_dPdF_incremental_fix;
#elif 0
        if(apply_definiteness_fix) rotated_dPdF.Make_Positive_Definite();
#endif

        return rotated_dPdF;
    }

};

}

#endif
