subroutine Jacobian(parameters,J)
    use nrtype; use dimensions; use optimization; use pdfs;use structural_p2;use grids;use targets
    implicit none
    real(SP),dimension(moment_conditions,parameters_to_est),intent(out)::J
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    real(SP),dimension(parameters_to_est,2)::parameters_new
    real(SP),dimension(moment_conditions,1)::moments,densities,moments_nom,densities_nom
    real(SP),dimension(moment_conditions,2)::moments_new,densities_new,moments_new_nom
    integer::i_l,p_l,m_l
    real(SP),dimension(parameters_to_est)::eps
    character::pause_k
    J=0.0_sp
    call SMM_se(parameters,moments,densities)
    do i_l=1,moment_conditions
        if (data_moments(i_l,1)==-9.0_sp) then
            moments(i_l,1)=-9.0_sp
            densities(i_l,1)=-9.0_sp
        end if
    end do
    
    call empty_missing(moments,moments_nom,int(moment_conditions),real_moments)
    call empty_missing(densities,densities_nom,int(moment_conditions),real_moments)
    
    eps=0.001_sp
    do p_l=1,parameters_to_est
        if (p_l==3) then
            eps(p_l)=0.01_sp
        end if
    end do
    
    do p_l=1,parameters_to_est
        do m_l=1,2
            parameters_new(:,m_l)=parameters
            if (m_l==1)then
                parameters_new(p_l,m_l)=parameters_new(p_l,m_l)+eps(p_l)*abs(parameters(p_l))
            else
                parameters_new(p_l,m_l)=parameters_new(p_l,m_l)-eps(p_l)*abs(parameters(p_l))
            end if
            call SMM_se(parameters_new(:,m_l),moments_new(:,m_l),densities_new(:,m_l))
            do i_l=1,moment_conditions
                if (data_moments(i_l,1)==-9.0_sp) then
                    moments_new(i_l,m_l)=-9.0_sp
                end if
            end do
            call empty_missing(moments_new(:,m_l),moments_new_nom(:,m_l),int(moment_conditions),real_moments)
        end do
        J(1:real_moments,p_l)=(moments_new_nom(1:real_moments,1)-moments_new_nom(1:real_moments,2))/(parameters_new(p_l,1)-parameters_new(p_l,2))*densities_nom(1:real_moments,1)
    end do        
      
end subroutine
    