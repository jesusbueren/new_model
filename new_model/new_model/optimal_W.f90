subroutine optimal_W()
    use simulation_input; use nrtype; use files_savings; use dimensions; use structural_p1; use structural_p2;use targets
    implicit none
    real(DP),dimension(moment_conditions,moment_conditions)::W_opt_d
    !Loop variables
    integer::i_l,t_l,m_l,ind,real_moments,s_l,i_l2,f_l
    !Individual data from HRS variables
    integer,dimension(indv,obs):: IC_q
    real(SP),dimension(indv,obs):: fc_h,NW,ic_h,beq100,govmd
    !Random draw if computing uncertainty of data
    real(SP)::u
    !Moment variables variables
    real(SP),dimension(L_PI,clusters)::fc_pi_h_it
    real(SP),dimension(L_PI,f_t)::beq100_ic_it
    real(SP),dimension(f_t,clusters)::fc_ic_h_it,md_ic_h_it
    real(SP),dimension(f_t,obs,groups)::nw_ic_it  
    real(SP),dimension(L_PI,obs,groups)::assets_pi_age_group_it,assets_pi_age_group_it_b
    real(SP),dimension(moment_conditions,obs)::moments_it
    real(SP),dimension(moment_conditions,indv)::moments_i,moments_i_new
    !Variance-covariance
    real(SP),dimension(moment_conditions,moment_conditions,samples_per_i)::Phi_s
    real(DP),dimension(moment_conditions,moment_conditions)::Phi_d
    real(SP),dimension(8,obs,indv)::data_moments_HRS
    !Charge individuals HRS data
    open(unit=9,file=fdir_inputs//'data_moments.txt')
        read(9,*) data_moments_HRS
    close(9)
    fc_h=reshape(data_moments_HRS(1,:,:), (/indv, obs/), order = (/ 2, 1 /))
    IC_q=reshape(data_moments_HRS(2,:,:), (/indv, obs/), order = (/ 2, 1 /))
    NW=reshape(data_moments_HRS(3,:,:), (/indv, obs/), order = (/ 2, 1 /))
    NW=NW/1000.0_sp
    beq100=reshape(data_moments_HRS(4,:,:), (/indv, obs/), order = (/ 2, 1 /))
    ic_h=reshape(data_moments_HRS(5,:,:), (/indv, obs/), order = (/ 2, 1 /))
    govmd=reshape(data_moments_HRS(6,:,:), (/indv, obs/), order = (/ 2, 1 /))
    
    !Compute data moments from original data
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    s_l=1
    Phi_s=0.0_sp
    do s_l=1,samples_per_i
        print*,s_l
        h_i=-9.0
        do i_l=1,indv;
            !Sample the family type 
            call RANDOM_NUMBER(u)
            if (u<IC_pr_i(i_l,1)) then
                f_l=1
            elseif (u<sum(IC_pr_i(i_l,1:2))) then
                f_l=2
            else
                f_l=3
            end if
            assets_pi_age_group_it=0.0_sp
            assets_pi_age_group_it_b=0.0_sp
            beq100_ic_it=0.0_sp
            nw_ic_it=0.0_sp
            fc_pi_h_it=0.0_sp
            fc_ic_h_it=0.0_sp
            md_ic_h_it=0.0_sp
            do t_l=obs,1,-1; 
                !Sample health status using Kim smoother
                if (s_h_i(i_l,t_l,1)/=-9.0_sp) then
                    call RANDOM_NUMBER(u)
                    ind=1
                    if (t_l==obs) then
                        do while (h_i(i_l,t_l,1)==-9 .and. ind<=clusters)
                            if (u<sum(s_h_i(i_l,t_l,1:ind)).or. ind==clusters) then
                                h_i(i_l,t_l,1)=ind
                            else
                                ind=ind+1
                            end if
                        end do
                    elseif (h_i(i_l,t_l+1,1)==-9) then
                        do while (h_i(i_l,t_l,1)==-9 .and. ind<=clusters)
                            if (u<sum(s_h_i(i_l,t_l,1:ind)).or. ind==clusters) then
                                h_i(i_l,t_l,1)=ind
                            else
                                ind=ind+1
                            end if
                        end do
                    else
                        do while (h_i(i_l,t_l,1)==-9 .and. ind<=clusters)
                            if (u<sum(f_h_i(i_l,t_l,1:ind)*H_av(1:ind,h_i(i_l,t_l+1,1),generation_i(i_l)+t_l-1,PI_q_i(i_l),gender_i(i_l)))/sum(H_av(1:clusters,h_i(i_l,t_l+1,1),generation_i(i_l)+t_l-1,PI_q_i(i_l),gender_i(i_l))*f_h_i(i_l,t_l,:)) .or. ind==clusters) then !Kim smoother
                                h_i(i_l,t_l,1)=ind
                            else
                                ind=ind+1
                            end if
                        end do
                    end if     
                else
                    h_i(i_l,t_l,1)=-9
                end if
            end do
            do t_l=1,obs
                if (fc_h(i_l,t_l)/=-9.0_sp) then
                    fc_pi_h_it(PI_q_i(i_l),h_i(i_l,t_l,1))=fc_h(i_l,t_l)-data_lfc_PI(PI_q_i(i_l),h_i(i_l,t_l,1))
                end if
                if (beq100(i_l,t_l)/=-9) then
                    beq100_ic_it(PI_q_i(i_l),f_l)=beq100(i_l,t_l)-data_beq100_ic(PI_q_i(i_l),f_l)
                end if
                if (NW(i_l,t_l)/=-9.0_sp) then
                    !By ic
                    if (NW(i_l,t_l)<data_NW_IC1(f_l,t_l,group_i(i_l))) then
                        nw_ic_it(f_l,t_l,group_i(i_l))=0.25_sp
                    else
                        nw_ic_it(f_l,t_l,group_i(i_l))=-0.75_sp
                    end if 
                end if
                    
                if (fc_h(i_l,t_l)/=-9.0_sp .and. IC_q(i_l,t_l)/=-9 .and. ic_h(i_l,t_l)/=-9.0_sp) then
                    fc_ic_h_it(IC_q(i_l,t_l),h_i(i_l,t_l,1))=fc_h(i_l,t_l)-data_lfc_IC(IC_q(i_l,t_l),h_i(i_l,t_l,1))
                end if
                if (govmd(i_l,t_l)/=-9.0_sp .and. IC_q(i_l,t_l)/=-9) then
                    md_ic_h_it(IC_q(i_l,t_l),h_i(i_l,t_l,1))=govmd(i_l,t_l)-data_govmd_IC(IC_q(i_l,t_l),h_i(i_l,t_l,1))
                end if
                if (NW(i_l,t_l)/=-9.0_sp ) then  !.and. (iwendy_i(i_l)==1998 .or. iwendy_i(i_l)==1999)
                    !By PIq
                    if (NW(i_l,t_l)<data_NW_PI1(PI_q_i(i_l),t_l,group_i(i_l))) then
                        assets_pi_age_group_it(PI_q_i(i_l),t_l,group_i(i_l))=0.5_sp
                    else
                        assets_pi_age_group_it(PI_q_i(i_l),t_l,group_i(i_l))=-0.5_sp
                    end if 
                end if

                moments_it(:,t_l)=(/reshape(assets_pi_age_group_it,(/L_PI*obs*groups,1/)),&
                                    reshape(beq100_ic_it,(/L_PI*f_t,1/)), &
                                    reshape(nw_ic_it,(/f_t*obs*groups,1/)), &
                                    reshape(fc_pi_h_it,(/L_PI*clusters,1/)), &
                                    reshape(fc_ic_h_it,(/f_t*clusters,1/)), &
                                    reshape(md_ic_h_it,(/f_t*clusters,1/))/)
            end do
            moments_i(:,i_l)=sum(moments_it,2)/real(obs)
            do m_l=1,moment_conditions
                if (data_moments(m_l,1)==-9.0_sp) then
                    moments_i(m_l,i_l)=-9.0_sp
                end if
            end do
            call empty_missing(moments_i(:,i_l),moments_i_new(:,i_l),int(moment_conditions),real_moments)
            Phi_s(1:real_moments,1:real_moments,s_l)=Phi_s(1:real_moments,1:real_moments,s_l)+matmul(moments_i_new(1:real_moments,i_l:i_l),transpose(moments_i_new(1:real_moments,i_l:i_l)))
        end do   
        Phi_s(1:real_moments,1:real_moments,s_l)=Phi_s(1:real_moments,1:real_moments,s_l)/real(indv)      
end do !simulation loop
Phi=sum(Phi_s,3)/real(samples_per_i)
Phi_d=dble(Phi)
call inverse(Phi_d(1:real_moments,1:real_moments),W_opt_d(1:real_moments,1:real_moments),real_moments)
W_opt=real(W_opt_d)
print*,'Diagonal of the optimal weigthing matrix'
do m_l=1,real_moments
    print*,m_l,W_opt(m_l,m_l)
end do


end subroutine