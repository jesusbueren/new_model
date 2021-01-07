subroutine compute_se(parameters)
    use dimensions; use nrtype; use optimization; use structural_p2; use simulation_input
    implicit none
    real(SP),dimension(parameters_to_est)::parameters
    real(SP),dimension(moment_conditions,parameters_to_est)::J
    real(DP),dimension(parameters_to_est,parameters_to_est)::V,inv_v
    integer::p_l,ns,i_l,j_l
    real(SP)::obj_fct_original
    real(SP),dimension(moment_conditions,1)::mean
    real(SP),dimension(moment_conditions,samples_per_i)::matrix_s_new
    interface
        function SMM(parameters)
        use nrtype;use dimensions
        implicit none
        real(SP),dimension(parameters_to_est),intent(in)::parameters
        real(SP)::SMM
        end function
    end interface
       
    !Compute Jacobian
    call Jacobian(parameters,J)
    
    !Variance Covariance matrix of estimated parameters
    V=matmul(matmul(transpose(dble(J(1:real_moments,1:parameters_to_est))),dble(W_opt(1:real_moments,1:real_moments))) &
             ,dble(J(1:real_moments,1:parameters_to_est)))
    call inverse(V,inv_V,parameters_to_est)
    !If using the optimal weighting matrix
    V=(1.0d0+1.0d0/dble(samples_per_i))/dble(indv)*inv_V
    
    print*,'standard errors'
    do p_l=1,parameters_to_est
        se(p_l)=real(sqrt(V(p_l,p_l)))
        print*,se(p_l)
    end do
    

end subroutine