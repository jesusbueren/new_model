subroutine charge_simulation_input_moments(data_NW_PI1,data_NW_PI,&
                                            data_NW_PI1b,data_NW_PIb,&
                                            data_beq100_IC,&
                                            data_NW_IC1,data_NW_IC,&
                                            data_lfc_PI,&
                                            data_lfc_IC,&
                                            data_NW_h_ut)
    use simulation_input; use nrtype; use files_savings; use dimensions; use structural_p1; use pdfs; use HRS_data
    implicit none
    real(SP),dimension(L_PI,f_t),intent(inout)::data_beq100_IC
    real(SP),dimension(L_PI,clusters),intent(inout)::data_lfc_PI
    real(SP),dimension(f_t,clusters),intent(inout)::data_lfc_IC
    real(SP),dimension(f_t,obs,groups),intent(inout)::data_NW_IC,data_NW_IC1
    real(SP),dimension(L_PI,obs,groups),intent(inout)::data_NW_PI,data_NW_PI1,data_NW_PI1b,data_NW_PIb
    real(SP),dimension(2,obs),intent(inout)::data_NW_h_ut
    real(SP),dimension(L_PI,obs,groups,samples_per_i)::data_NW_PI1_ns,data_NW_PI1_ns_b
    real(SP),dimension(2,obs,samples_per_i)::data_NW_h_ut_ns
    !Loop variables
    integer::i_l,t_l,pi_l,h_l,g_l,f_l,s_l,ind,ic_p,ic_l,q_l
    !Counter variables
    integer,dimension(2,obs)::counter_ut
    integer,dimension(L_PI,clusters)::counter_pi_h
    integer,dimension(L_PI,clusters)::counter_pi_h2
    integer,dimension(f_t,clusters)::counter_ic_h
    integer,dimension(f_t,obs,groups)::counter_ic_nw
    integer,dimension(L_PI,f_t)::counter_beq100
    integer,dimension(L_PI,obs,groups)::counter_pi_age_group,counter_pi_age_group_b
    !hours of care variables
    real(SP),dimension(L_PI,clusters,6500)::fc_pi_h
    real(SP),dimension(f_t,clusters,15000)::fc_ic_h,ic_ic_h
    integer,dimension(L_PI,f_t,10000)::beq100_ic
    !Store vector variables    
    real(SP),dimension(L_PI,obs,groups,1500)::assets_pi_age_group,assets_pi_age_group_b !1500: maximum number of individuals in a group
    real(SP),dimension(2,obs,8000)::assets_ut
    real(SP),dimension(f_t,obs,groups,4000)::assets_ic    
    !Across latent health types
    real(SP),dimension(L_PI,clusters,samples_per_i)::data_NW_PI_s
    real(SP),dimension(L_PI,clusters,samples_per_i)::data_lfc_PI_s
    real(SP),dimension(f_t,clusters,samples_per_i)::data_lfc_IC_s,l_ic_s
    real(SP),dimension(L_PI,clusters,samples_per_i)::pr_pi_h_s,pr_pi_h2_s
    real(SP),dimension(f_t,clusters,samples_per_i)::pr_ic_h_s,pr_ic_h2_s
    real(SP),dimension(f_t,obs,groups,samples_per_i)::assets_ic_ns
    real(SP),dimension(L_PI,f_t,samples_per_i)::data_beq100_IC_s
    !Random draw
    real(SP)::u
    !Moments
    real(SP),dimension(5,obs,indv)::data_moments
    
    !Charge individuals HRS data
    open(unit=9,file=fdir_inputs//'data_moments.txt')
        read(9,*) data_moments
    close(9)
    fc_h=reshape(data_moments(1,:,:), (/indv, obs/), order = (/ 2, 1 /))
    IC_q=reshape(data_moments(2,:,:), (/indv, obs/), order = (/ 2, 1 /))
    NW=reshape(data_moments(3,:,:), (/indv, obs/), order = (/ 2, 1 /))
    NW=NW/1000.0_sp
    beq100=reshape(data_moments(4,:,:), (/indv, obs/), order = (/ 2, 1 /))
    ic_h=reshape(data_moments(5,:,:), (/indv, obs/), order = (/ 2, 1 /))
    
    !Compute data moments from original data
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    do s_l=1,samples_per_i
        counter_pi_h=0
        counter_pi_h2=0
        counter_ic_h=0
        counter_pi_age_group=0
        counter_pi_age_group_b=0
        counter_ut=0
        fc_pi_h=0.0_sp
        fc_ic_h=0.0_sp
        ic_ic_h=0.0_sp
        beq100_ic=0
        counter_ic_nw=0
        counter_beq100=0
        h_i=-9.0
        !Store vector with all: 
        !       - Formal care hours across PI quartiles and health status
        !       - Medicaid recipients across PI quartiles and health status
        !       - Formal care and informal care hours across observed families and health status
        !       - Wealth across PI quartiles and age for individuals observed in 1998 or 1999
        do i_l=1,indv;
            !Sample the family type 
            call RANDOM_NUMBER(u)
            if (u<IC_pr_i(i_l,1)) then
                f_l=1
            elseif (u<sum(IC_pr_i(i_l,1:2))) then
                f_l=2
            else
                print*,'error in charge sim input'
            end if
            do t_l=obs,1,-1
                !Sample health status using Kim Smoother
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
                    elseif ( h_i(i_l,t_l+1,1)==-9) then
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
                    counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,1))=counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,1))+1
                    fc_pi_h(PI_q_i(i_l),h_i(i_l,t_l,1),counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,1)))=fc_h(i_l,t_l)
                end if
                if (beq100(i_l,t_l)/=-9) then
                    counter_beq100(PI_q_i(i_l),f_l)=counter_beq100(PI_q_i(i_l),f_l)+1
                    beq100_ic(PI_q_i(i_l),f_l,counter_beq100(PI_q_i(i_l),f_l))=beq100(i_l,t_l)
                end if
                if (fc_h(i_l,t_l)/=-9.0_sp .and. IC_q(i_l,t_l)/=-9 .and. ic_h(i_l,t_l)/=-9.0_sp) then
                    counter_ic_h(IC_q(i_l,t_l),h_i(i_l,t_l,1))=counter_ic_h(IC_q(i_l,t_l),h_i(i_l,t_l,1))+1
                    fc_ic_h(IC_q(i_l,t_l),h_i(i_l,t_l,1),counter_ic_h(IC_q(i_l,t_l),h_i(i_l,t_l,1)))=fc_h(i_l,t_l)
                    ic_ic_h(IC_q(i_l,t_l),h_i(i_l,t_l,1),counter_ic_h(IC_q(i_l,t_l),h_i(i_l,t_l,1)))=ic_h(i_l,t_l)
                end if
                if (NW(i_l,t_l)/=-9.0_sp ) then
                    counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l))=counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l))+1
                    assets_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l),counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l)))=NW(i_l,t_l)
                    if (h_i(i_l,t_l,1)==1) then
                        counter_ut(1,t_l)=counter_ut(1,t_l)+1
                        assets_ut(1,t_l,counter_ut(1,t_l))=NW(i_l,t_l)
                    else
                        counter_ut(2,t_l)=counter_ut(2,t_l)+1
                        assets_ut(2,t_l,counter_ut(2,t_l))=NW(i_l,t_l)
                    end if
                end if
                if (NW(i_l,t_l)/=-9.0_sp) then
                    counter_ic_nw(f_l,t_l,group_i(i_l))=counter_ic_nw(f_l,t_l,group_i(i_l))+1
                    assets_ic(f_l,t_l,group_i(i_l),counter_ic_nw(f_l,t_l,group_i(i_l)))=NW(i_l,t_l)
                end if
            end do
        end do
    
        !Compute: 
        !       - Average formal care hours across PI quartiles and health status
        !       - Medicaid recipiency rates across PI quartiles and health status
        !       - Average formal care and informal care hours across observed families and health status
        !       - Median wealth across PI quartiles and age for individuals observed in 1998 or 1999
    
        data_lfc_PI_s(:,:,s_l)=sum(fc_pi_h,3)/counter_pi_h
        data_lfc_PI_s(:,1,s_l)=-9.0_sp
        data_beq100_IC_s(:,:,s_l)=real(sum(beq100_ic,3))/counter_beq100
        data_lfc_IC_s(:,:,s_l)=sum(fc_ic_h,3)/counter_ic_h
        data_lfc_IC_s(:,1,s_l)=-9.0_sp
    
        !Compute hours of informal care across families l_ic_s(3,:,s_l)
        l_ic_s(:,:,s_l)=sum(ic_ic_h,3)/counter_ic_h
        l_ic_s(:,:,s_l)=l_ic_s(:,:,s_l)*2.0_sp*365.0_sp
        l_ic_s(:,1,s_l)=0.0_sp
        
        !Wealth moments by PIq
        do pi_l=1,L_PI; do t_l=1,7; do g_l=1,4
            if (counter_pi_age_group(pi_l,t_l,g_l)>min_obs) then
                call compute_percentile(assets_pi_age_group(pi_l,t_l,g_l,1:counter_pi_age_group(pi_l,t_l,g_l)), &
                                    counter_pi_age_group(pi_l,t_l,g_l),50, &
                                    data_NW_PI1_ns(pi_l,t_l,g_l,s_l)) 
            end if 
            if (counter_pi_age_group_b(pi_l,t_l,g_l)>min_obs) then
                call compute_percentile(assets_pi_age_group_b(pi_l,t_l,g_l,1:counter_pi_age_group_b(pi_l,t_l,g_l)), &
                                    counter_pi_age_group_b(pi_l,t_l,g_l),50, &
                                    data_NW_PI1_ns_b(pi_l,t_l,g_l,s_l)) 
            end if 
        end do; end do; end do
        !Wealth moments by h 
        do h_l=1,2; do t_l=1,7
            if (counter_ut(h_l,t_l)>min_obs) then
                call compute_percentile(assets_ut(h_l,t_l,1:counter_ut(h_l,t_l)), &
                                    counter_ut(h_l,t_l),50, &
                                    data_NW_h_ut_ns(h_l,t_l,s_l)) 
            end if 
        end do; end do
        !Wealth moments by ic
        do f_l=1,f_t; do t_l=1,7; do g_l=1,4
            if (counter_ic_nw(f_l,t_l,g_l)>min_obs) then
                call compute_percentile(assets_ic(f_l,t_l,g_l,1:counter_ic_nw(f_l,t_l,g_l)), &
                                    counter_ic_nw(f_l,t_l,g_l),75, &
                                    assets_ic_ns(f_l,t_l,g_l,s_l)) 
            end if 
        end do;end do;end do
        
    end do
    data_lfc_PI=sum(data_lfc_PI_s,3)/real(samples_per_i)
    data_beq100_IC=sum(data_beq100_IC_s,3)/real(samples_per_i)
    l_ic=sum(l_ic_s,3)/real(samples_per_i) !l_ic
    data_lfc_IC=sum(data_lfc_IC_s,3)/real(samples_per_i)

    !Wealth moments by PIq
    do pi_l=1,L_PI; do t_l=1,7; do g_l=1,4
        if (counter_pi_age_group(pi_l,t_l,g_l)>min_obs) then
            data_NW_PI1(pi_l,t_l,g_l)=sum(data_NW_PI1_ns(pi_l,t_l,g_l,:))/real(samples_per_i)
            data_NW_PI(pi_l,t_l,g_l)=0.5_sp
        end if  
        if (counter_pi_age_group_b(pi_l,t_l,g_l)>min_obs) then
            data_NW_PI1b(pi_l,t_l,g_l)=sum(data_NW_PI1_ns_b(pi_l,t_l,g_l,:))/real(samples_per_i)
            data_NW_PIb(pi_l,t_l,g_l)=0.5_sp
        end if  
    end do; end do; end do
    !Wealth moments by h
    do h_l=1,2; do t_l=1,7;
        if (counter_ut(h_l,t_l)>min_obs) then
            data_NW_h_ut(h_l,t_l)=sum(data_NW_h_ut_ns(h_l,t_l,:))/real(samples_per_i)
        end if  
    end do; end do
    !Wealth moments by family
    do f_l=1,f_t;do t_l=1,7; do g_l=1,4
        if (counter_ic_nw(f_l,t_l,g_l)>min_obs) then
            data_NW_IC1(f_l,t_l,g_l)=sum(assets_ic_ns(f_l,t_l,g_l,:))/real(samples_per_i)
            data_NW_IC(f_l,t_l,g_l)=0.75_sp
        end if  
    end do;end do;end do

end subroutine