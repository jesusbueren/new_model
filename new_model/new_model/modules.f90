MODULE nrtype
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.true.)
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
	REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
	REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
	REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
	TYPE sprs2_sp
		INTEGER(I4B) :: n,len
		REAL(SP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_sp
	TYPE sprs2_dp
		INTEGER(I4B) :: n,len
		REAL(DP), DIMENSION(:), POINTER :: val
		INTEGER(I4B), DIMENSION(:), POINTER :: irow
		INTEGER(I4B), DIMENSION(:), POINTER :: jcol
	END TYPE sprs2_dp
END MODULE nrtype
    
module dimensions
use nrtype
    implicit none
    !nkk: grid point cash on hand (coh) grid
    !f_t:family type
    !h_t:health type
    !g_t: gender type
    !PI_t: permanent income type                                                     
    integer,parameter::nkk=300,f_t=2,clusters=4,L_gender=2,L_PI=5,generations=21,nzz=3,nzz2=1,variables=12,groups=4,parameters_to_est=8,obs=9, &
                        wealth_q=3,L_PI2=10,samples_per_i=50,min_obs=39
    integer,parameter:: moment_conditions=L_PI*obs*groups*2+L_PI*clusters+f_t*obs*groups+L_PI*clusters+f_t*clusters
end module dimensions
        
module grids
    use dimensions; use nrtype
    implicit none
    real(SP)::coh_min=0.0_sp,coh_max=600.0_sp
    real(SP),dimension(nkk)::coh_grid 
    real(SP),dimension(nzz,1)::ep_grid 
    real(SP),dimension(nzz,1)::xi_grid 
    real(SP),dimension(nzz,1)::pr_ep,pr0_xi
    real(SP),dimension(nzz,nzz)::pr_pxi
end module grids
    
module targets
    use dimensions; use nrtype
    implicit none
    real(SP),dimension(L_PI,obs,groups)::data_NW_PI,data_NW_PI1,data_NW_PIb,data_NW_PI1b
    real(SP),dimension(2,obs)::data_NW_h_ut
    real(SP),dimension(L_PI,clusters)::data_MD_PI
    real(SP),dimension(L_PI,clusters)::data_lfc_PI
    real(SP),dimension(f_t,clusters)::data_lfc_IC
    real(SP),dimension(f_t,obs,groups)::data_NW_IC,data_NW_IC1
    real(SP),dimension(moment_conditions,1)::data_moments_new,data_moments1,data_moments
end module targets
    
module structural_p2
    use dimensions; use nrtype
    implicit none
    real(SP),dimension(nzz2,1)::pr_varep
    real(SP),dimension(nzz2,clusters-1)::varep_grid
    real(SP):: sigma,sigma_beq,nu,beta,omega
    real(SP),dimension(clusters)::sigma_varep
    real(SP),dimension(clusters):: kappa_h,delta_h,u_bar_no_f,x_bar
    real(SP),dimension(clusters,f_t,nzz2)::u_bar
    real(SP),dimension(clusters-1)::alpha_mu,sigma2_varep
    real(SP),dimension(clusters,nzz2)::c_bar,l_bar
    real(SP),dimension(f_t):: lambda,delta
    real(SP),dimension(clusters,nzz2):: mu
    real(SP),dimension(moment_conditions,moment_conditions):: W_opt,Phi
    real(SP),dimension(parameters_to_est)::se
end module structural_p2
    
module structural_p1
    use dimensions; use nrtype
    implicit none
    integer,parameter::DIM=12+11+3
    real(SP),dimension(DIM,1)::med_coef
    real(SP)::rho,sigma2_ep,sigma2_ze
    real(SP),parameter::r=1.02_sp**2_sp-1.0_sp
    real(SP)::p_fc=18.0_sp/1000.0_sp
    real(SP),parameter::p_or=18.0_sp/1000.0_sp
    real(SP),dimension(L_PI2,L_gender)::b
    real(SP),dimension(f_t,clusters)::l_ic,l_ic_or
    real(SP),dimension(generations,L_gender,L_PI2,clusters,nzz,nzz)::m_exp_all,m_exp_all_or
    real(SP),dimension(clusters+1,clusters+1,generations,L_PI2,L_gender)::H_av
    real(SP)::load_ltci=0.32_sp,benefit_LTCI=5.0_sp !hours of care per day provided by insurance
    real(SP),dimension(L_PI,2,clusters)::price_ltci
end module structural_p1
    
module files_savings
implicit none
    character(LEN=76)::fdir_med="C:\Users\jbueren\Google Drive\JMP\Code\Health Dynamics\medical_expenditures\"
    character(LEN=64)::fdir_inputs="C:\Users\jbueren\Google Drive\JMP\Code\Structural Model\moments\"
    character(LEN=63)::path_h="C:\Users\jbueren\Google Drive\JMP\Code\Health Dynamics\Results\"
    !character(LEN=74)::fdir_med="C:\Users\Jesus\Google Drive\JMP\Code\Health Dynamics\medical_expenditures\"
    !character(LEN=62)::fdir_inputs="C:\Users\Jesus\Google Drive\JMP\Code\Structural Model\moments\"
    end module
    
module simulation_input
use nrtype; use dimensions
implicit none
    integer,parameter::indv=8472
    integer,dimension(indv):: gender_i, generation_i, PI_q_i, PI_q_i2, group_i,iwendy_i
    real(SP),dimension(indv,f_t):: IC_pr_i
    real(SP),dimension(indv):: x_i
    real(SP),dimension(indv,obs,clusters):: s_h_i,f_h_i
    real(SP),dimension(indv,obs):: dead_i
    integer,dimension(indv,obs,2):: h_i
end module
    
module optimization
use nrtype; use structural_p2; use dimensions
implicit none
    real(SP)::obj_fct_s
    real(SP),dimension(moment_conditions,1)::r1,f
    integer::real_moments
    real(SP),dimension(moment_conditions,moment_conditions)::matrix_s
    real(SP),dimension(L_PI,obs,groups)::moments_NW_PI1
    real(SP),dimension(f_t,obs,groups)::moments_NW_IC1
    real(SP),dimension(wealth_q,obs,groups)::moments_NW_all1
    real(SP):: obj_fct_min=1.0_sp/0.0_sp
end module
    
module pdfs
    use nrtype;use dimensions;use simulation_input
    implicit none
    real(SP),dimension(L_PI,obs,groups,nkk)::pdf_all_NW_PI
end module
    
module MD_reform
use dimensions;use nrtype
implicit none
real(SP),dimension(nkk,clusters+1,nzz,nzz2,L_gender,L_PI2,f_t)::V_70_or,V_70_new,V_70
integer::ind_or=1
real(SP)::p_sub=0.0_sp
end module

    
