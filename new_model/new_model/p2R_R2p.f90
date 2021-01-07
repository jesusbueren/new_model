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
                 parameters(8)/) !bequest intensity
                 !log(parameters(14)),& !scale 
                 !log(parameters(15)),& !scale 
                 !log(parameters(16)),& !scale
                 !log(parameters(17)),& !scale
                 !log(parameters(14)/(1.0_sp-parameters(14))), &!beta: anti logit function
                 !log(parameters(11)/(1.0_sp-parameters(13))), & !omega: anti logit function   
                 !log(parameters(17)), &
                 !log(parameters(18)), &
                 !log(parameters(19))/) 
    
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
    !parameters_n(14:17)=exp(parameters(14:17))
    !parameters_n(15)=exp(parameters(15))/(exp(parameters(15))+1.0_sp)
    !parameters_n(16:19)=exp(parameters(16:19))
    !parameters_n(13)=exp(parameters(13))/(exp(parameters(13))+1.0_sp)
    
    
    !Set parameters
    sigma=parameters_n(1)
    x_bar(1)=parameters_n(2)
    alpha_mu(1:3)=parameters_n(3:5)
    delta(1:f_t)=parameters_n(6)
    nu=parameters_n(7)
    lambda(1:f_t)=exp(parameters_n(8))

    sigma_beq=sigma
    beta=0.95_sp 
    omega=1.0_sp
    sigma_varep(1:4)=0.0000001
    sigma2_varep(1:3)=0.0_sp
    kappa_h=0.0_sp
    delta_h=1.0_sp
    
    !Utility at the floor when healthy
    do h_l=1,1
        if (h_l==1) then
            mu_av=0.0_sp
        elseif (h_l==2) then
            mu_av=exp(alpha_mu(1)+sigma2_varep(1)/2.0_sp)
        elseif (h_l==3) then
            mu_av=exp(alpha_mu(2)+sigma2_varep(2)/2.0_sp)
        elseif (h_l==4) then
            mu_av=exp(alpha_mu(3)+sigma2_varep(3)/2.0_sp)
        end if
        call solve_intratemporal_av(p_fc,x_bar(h_l),h_l,mu_av,0.0_sp,u_bar_no_f(h_l),l_fc,c)
    end do
    u_bar_no_f=u_bar_no_f(1)
    
end subroutine    
    