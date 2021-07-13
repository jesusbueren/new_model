subroutine p2R(parameters)
    use nrtype;use dimensions
    implicit none
    real(SP),dimension(parameters_to_est),intent(inout)::parameters
    
    parameters=(/log(parameters(1)-1.0_sp), & !RRA consumption
                 log(parameters(2)), & !av money transfer at floor h=1
                 parameters(3), & !shifter h=2
                 parameters(4), & !shifter h=3
                 parameters(5), & !shifter h=4
                 log(parameters(6)), & !beq curvature 
                 log(parameters(7)-1.0_sp),& !RRA LTC 
                 parameters(8), & !bequest intensity (F=1)
                 log(parameters(9)/(1.0_sp-parameters(9))), & !share from CES
                 log(1.0_sp-parameters(10))/) !subtitution parameter

    
end subroutine
    
subroutine R2p(parameters,parameters_n)
    use nrtype;use dimensions; use structural_p2; use structural_p1
    implicit none
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    real(SP),dimension(parameters_to_est),intent(out)::parameters_n
    integer::h_l,f_l
    real(SP)::mu_av,l_fc,c
    
    !Tranform parameters from the real line
    parameters_n(1)=exp(parameters(1))+1.0_sp !RRA consumption
    parameters_n(2)=exp(parameters(2)) !av money transfer at floor h=1
    parameters_n(3:5)=parameters(3:5) !shifter h=2,3,4
    parameters_n(6)=exp(parameters(6))  !beq curvature 
    parameters_n(7)=exp(parameters(7))+1.0_sp !RRA LTC 
    parameters_n(8)=parameters(8)  !bequest intensity
    parameters_n(9)=1.0_sp/(1.0_sp+exp(-parameters(9))) !share parameter
    parameters_n(10)=1.0_sp-exp(parameters(10))!share parameter
    
    !Set parameters
    sigma=parameters_n(1)
    x_bar(1)=parameters_n(2)
    alpha_mu(1:3)=parameters_n(3:5)
    delta(1:f_t)=parameters_n(6)
    nu=parameters_n(7)
    lambda(1:f_t)=exp(parameters_n(8))
    share_p=parameters_n(9)
    subs_p=parameters_n(10)

    sigma_beq=sigma
    beta=0.95_sp 
    omega=1.0_sp
    kappa_h=0.0_sp
    delta_h=1.0_sp
    
    !Utility at the floor when healthy
    do h_l=1,1
        call solve_intratemporal(p_fc(1),x_bar(h_l),h_l,0.0_sp,u_bar_no_f(h_l),l_fc,c)
    end do
    u_bar_no_f=u_bar_no_f(1)
    
end subroutine    
    