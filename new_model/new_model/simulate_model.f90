subroutine simulate_model(a_policy,g_policy,lfc_x,beq100_policy, &
                          model_moments1,model_moments)
    use dimensions;use grids; use nrtype; use simulation_input; use structural_p1; use pdfs; use structural_p2;use targets; use HRS_data
    implicit none
    real(SP),dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations),intent(in)::beq100_policy
    integer,dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations),intent(in)::a_policy,g_policy
    real(SP),dimension(nkk,clusters,f_t,L_PI2),intent(in)::lfc_x
    real(SP),dimension(moment_conditions,1),intent(out)::model_moments1,model_moments
    real(SP),dimension(2,obs)::model_NW_h_ut
    !Loop variables
    integer::i_l,t_l,pos_x,xi_l,xi_l2,ts_l2,ind,k2_l,ns,g_l,pi_l,ic_l,h_l,f_l2,nwq_l,s_l,it,f_ll
    integer,dimension(2)::f_l
    !Counter variables
    integer,dimension(L_PI,obs,groups)::counter_pi_age_group,counter_pi_age_group_b
    integer,dimension(f_t,2)::counter_ic_age_group,counter_ic_age_group_ut
    integer,dimension(obs,groups)::counter_all_age_group
    integer,dimension(L_PI,clusters)::counter_pi_h
    integer,dimension(L_PI,clusters)::counter_pi_h2
    integer,dimension(f_t,clusters)::counter_ic_h
    integer,dimension(f_t)::counter_ic
    integer,dimension(2,obs)::counter_ut
    integer,dimension(f_t,obs,groups)::counter_ic_nw
    integer,dimension(L_PI,f_t)::counter_beq100
    !Hours of care variables and medicaid variables
    real(SP),dimension(L_PI,clusters,indv)::lfc_pi_h
    real(SP),dimension(L_PI,clusters,indv)::MD_pi_h
    real(SP),dimension(f_t,clusters,indv*2)::lfc_ic_h,MD_ic_h
    real(SP),dimension(L_PI,f_t,samples_per_i)::av_beq100_ic
    real(SP),dimension(L_PI,clusters,samples_per_i)::av_lfc_pi_h
    real(SP),dimension(f_t,clusters,samples_per_i)::av_lfc_ic_h
    real(SP),dimension(2,obs,samples_per_i)::model_NW_h_ut_ns
    real(SP),dimension(f_t,obs,groups,samples_per_i)::assets_ic_ns
    !Store vector variables for assets
    real(SP),dimension(L_PI,obs,groups,1400)::assets_pi_age_group,assets_pi_age_group_b !1400: maximum number of individuals in a group
    real(SP),dimension(2,obs,6000)::assets_ut
    real(SP),dimension(L_PI,f_t,30000)::beq100_ic
    real(SP),dimension(f_t,obs,groups,4000)::assets_ic
    !Model moments
    real(SP),dimension(L_PI,obs,groups)::moments_NW_PI1,moments_NW_PI1b
    real(SP),dimension(L_PI,f_t)::moments_beq100_IC
    real(SP),dimension(L_PI,clusters)::moments_MD_PI
    real(SP),dimension(L_PI,clusters)::moments_lfc_PI
    real(SP),dimension(f_t,clusters)::moments_lfc_IC
    real(SP),dimension(f_t,obs,groups)::moments_NW_IC1
    real(SP),dimension(L_PI,obs,groups,samples_per_i)::med_assets_pi_age_group,med_assets_pi_age_group_b
    !Moment variables variables
    real(SP),dimension(L_PI,clusters)::fc_pi_h_it
    real(SP),dimension(L_PI,clusters)::govmd_pi_h_it
    real(SP),dimension(f_t,obs,groups)::nw_ic_it
    real(SP),dimension(f_t,clusters)::fc_ic_h_it
    real(SP),dimension(L_PI,f_t)::beq100_it
    real(SP),dimension(L_PI,obs,groups)::assets_pi_age_group_it,assets_pi_age_group_it_b
    real(SP),dimension(moment_conditions,obs,indv)::moments_it
    real(SP),dimension(moment_conditions,samples_per_i)::moments_s
    !Simulation variables
    real(SP)::u,k2,x_pr_md
    real(SP),dimension(obs+1)::x_it,g_it,a_it
    integer,dimension(1)::seed=254
    character::pause_k
    !Gvt transfers
    real(SP),dimension(L_PI)::gvt_pi
    !Set the seed
    call random_seed(PUT=seed) 
    !Initialize value for moments
    moments_s=0.0_sp
    med_assets_pi_age_group=-9.0_sp
    med_assets_pi_age_group_b=-9.0_sp
    av_beq100_ic=-9.0_sp
    av_lfc_pi_h=-9_sp
    av_lfc_ic_h=-9_sp
    do ns=1,samples_per_i
        !Initialize value for moments in each simulation
        counter_pi_age_group=0
        counter_pi_age_group_b=0
        counter_ic_h=0
        counter_ic=0
        counter_pi_h=0
        counter_pi_h2=0
        counter_ut=0
        counter_ic_nw=0
        counter_beq100=0
        assets_pi_age_group=-9.0_sp
        assets_pi_age_group_b=-9.0_sp
        assets_ic=-9.0_sp
        lfc_pi_h=-9.0_sp
        lfc_ic_h=0.0_sp
        moments_it=0.0_sp
        MD_pi_h=-9
        MD_ic_h=-9
        do i_l=1,indv
            a_it(1)=x_i(i_l)
            !Sample the true and measured family type 
            do it=1,2
                call RANDOM_NUMBER(u)
                if (u<IC_pr_i(i_l,1)) then
                    f_l(it)=1
                elseif (u<sum(IC_pr_i(i_l,1:2))) then
                    f_l(it)=2
                else
                    print*,'error in sim model'
                end if
            end do
            !Sample true health and measured health using Kim smoother
            do s_l=1,2; do t_l=obs,1,-1
                if (s_h_i(i_l,t_l,1)/=-9.0_sp) then
                    call RANDOM_NUMBER(u)
                    ind=1
                    if (t_l==obs) then
                        do while (h_i(i_l,t_l,s_l)==-9 .and. ind<=clusters)
                            if (u<sum(s_h_i(i_l,t_l,1:ind)).or. ind==clusters) then
                                h_i(i_l,t_l,s_l)=ind
                            else
                                ind=ind+1
                            end if
                        end do
                    elseif ( h_i(i_l,t_l+1,s_l)==-9) then
                        do while (h_i(i_l,t_l,s_l)==-9 .and. ind<=clusters)
                            if (u<sum(s_h_i(i_l,t_l,1:ind)).or. ind==clusters) then
                                h_i(i_l,t_l,s_l)=ind
                            else
                                ind=ind+1
                            end if
                        end do
                    else 
                        do while (h_i(i_l,t_l,s_l)==-9 .and. ind<=clusters)
                            if (u<sum(f_h_i(i_l,t_l,1:ind)*H_av(1:ind,h_i(i_l,t_l+1,s_l),generation_i(i_l)+t_l-1,PI_q_i2(i_l),gender_i(i_l)))/sum(H_av(1:clusters,h_i(i_l,t_l+1,s_l),generation_i(i_l)+t_l-1,PI_q_i2(i_l),gender_i(i_l))*f_h_i(i_l,t_l,:)) .or. ind==clusters) then !Kim smoother
                                h_i(i_l,t_l,s_l)=ind
                            else
                                ind=ind+1
                            end if
                        end do
                    end if     
                else
                    h_i(i_l,t_l,s_l)=-9
                end if
            end do;end do
            !Simulate decisions
            assets_pi_age_group_it=0.0_sp
            assets_pi_age_group_it_b=0.0_sp
            govmd_pi_h_it=0.0_sp
            nw_ic_it=0.0_sp
            fc_pi_h_it=0.0_sp
            fc_ic_h_it=0.0_sp
            do t_l=1,obs
                if (h_i(i_l,t_l,1) /= -9) then
                    !Sample initial medical persistent shock to characterize the state variable
                    call RANDOM_NUMBER(u)
                    if (t_l==1) then
                        xi_l=-9
                        ind=1
                        do while (xi_l==-9 .and. ind<=nzz)
                            if (u<sum(pr0_xi(1:ind,1)).or. ind==nzz) then
                                xi_l=ind
                            else
                                ind=ind+1
                            end if
                        end do
                            if (xi_l==-9) then
                                print*,'something wrong 0'
                                xi_l=nzz
                            end if
                    end if     
                    
                    !Sample persistent & transitory medical shock to simulate expenditures
                    xi_l2=-9
                    ind=1
                    call RANDOM_NUMBER(u)
                    do while (xi_l2==-9 .and. ind<=nzz)
                        if (u<sum(pr_pxi(xi_l,1:ind)) .or. ind==nzz) then
                            xi_l2=ind
                            xi_l=xi_l2
                        else
                            ind=ind+1
                        end if
                    end do
                    if (xi_l2==-9) then
                        print*,'something wrong 1'
                        xi_l2=nzz
                        xi_l=xi_l2
                    end if
                    ts_l2=-9
                    ind=1
                    call RANDOM_NUMBER(u)
                    do while (ts_l2==-9 .and. ind<=nzz)
                        if (u<sum(pr_ep(1:ind,1)) .or. ind==nzz) then
                            ts_l2=ind
                        else
                            ind=ind+1
                        end if
                    end do
                    if (ts_l2==-9) then
                        print*,'something wrong 2'
                        ts_l2=nzz
                    end if
                    
                    !Cash on hand
                    x_it(t_l)=(1+r)*a_it(t_l)-m_exp_all(generation_i(i_l)+t_l-1,gender_i(i_l),PI_q_i2(i_l),h_i(i_l,t_l,1),xi_l2,ts_l2)+b(PI_q_i2(i_l),gender_i(i_l))
                    pos_x=int(max(min(x_it(t_l),coh_grid(nkk)),0.0_sp)/(coh_grid(2)-coh_grid(1))+1.00000001_sp)
                    !Draw on the discrete choice
                    call RANDOM_NUMBER(u)

                    if (u<g_policy(pos_x,h_i(i_l,t_l,1),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l(1),generation_i(i_l)+t_l-1)) then
                        g_it(t_l)=1.0_sp
                        k2_l=1
                    else
                        g_it(t_l)=0.0_sp
                        k2_l=a_policy(pos_x,h_i(i_l,t_l,1),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l(1),generation_i(i_l)+t_l-1)
                    end if
                    if (t_l<obs) then
                        a_it(t_l+1)=coh_grid(k2_l)
                    end if

                    !Medicaid and hours of care
                
                    counter_pi_h2(PI_q_i(i_l),h_i(i_l,t_l,2))=counter_pi_h2(PI_q_i(i_l),h_i(i_l,t_l,2))+1
                    
                    counter_ic(f_l(2))=counter_ic(f_l(2))+1
                    if (beq100(i_l,t_l)/=-9.0_sp) then
                        counter_beq100(PI_q_i(i_l),f_l(2))=counter_beq100(PI_q_i(i_l),f_l(2))+1
                        beq100_ic(PI_q_i(i_l),f_l(2),counter_beq100(PI_q_i(i_l),f_l(2)))=beq100_policy(pos_x,h_i(i_l,t_l,1),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l(1),generation_i(i_l)+t_l-1)
                        beq100_it(PI_q_i(i_l),f_l(2))=beq100_ic(PI_q_i(i_l),f_l(2),counter_beq100(PI_q_i(i_l),f_l(2)))-data_beq100_IC(PI_q_i(i_l),f_l(2))
                    end if
                    
                    !Moments
                    !Store moments
                   ! if (iwendy_i(i_l)==1998 .or. iwendy_i(i_l)==1999) then  !if (iwendy_i(i_l)>1999) then
                    if (NW(i_l,t_l)/=-9.0_sp ) then
                        counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l))=counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l))+1
                        assets_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l),counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l)))=a_it(t_l)
                        if (a_it(t_l)<=data_NW_PI1(PI_q_i(i_l),t_l,group_i(i_l))) then
                            assets_pi_age_group_it(PI_q_i(i_l),t_l,group_i(i_l))=0.5_sp
                        else
                            assets_pi_age_group_it(PI_q_i(i_l),t_l,group_i(i_l))=-0.5_sp
                        end if 
                        counter_ic_nw(f_l(2),t_l,group_i(i_l))=counter_ic_nw(f_l(2),t_l,group_i(i_l))+1
                        assets_ic(f_l(2),t_l,group_i(i_l),counter_ic_nw(f_l(2),t_l,group_i(i_l)))=a_it(t_l)
                        if (a_it(t_l)<=data_NW_IC1(f_l(2),t_l,group_i(i_l))) then
                            nw_ic_it(f_l(2),t_l,group_i(i_l))=0.25_sp
                        else
                            nw_ic_it(f_l(2),t_l,group_i(i_l))=-0.75_sp
                        end if  
                    end if
                                                          
                    if (fc_h(i_l,t_l)/=-9.0_sp .and. IC_q(i_l,t_l)/=-9 .and. ic_h(i_l,t_l)/=-9.0_sp) then
                        counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2))=counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2))+1
                        counter_ic_h(f_l(1),h_i(i_l,t_l,2))=counter_ic_h(f_l(1),h_i(i_l,t_l,2))+1
                        if (g_it(t_l)==0) then
                        lfc_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2),counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2)))=lfc_x(pos_x-k2_l+1,h_i(i_l,t_l,1),f_l(1),PI_q_i2(i_l))/2.0_sp/365.0_sp
                        lfc_ic_h(f_l(1),h_i(i_l,t_l,2),counter_ic_h(f_l(1),h_i(i_l,t_l,2)))=          lfc_x(pos_x-k2_l+1,h_i(i_l,t_l,1),f_l(1),PI_q_i2(i_l))/2.0_sp/365.0_sp
                        else
                            lfc_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2),counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2)))=l_bar(h_i(i_l,t_l,1))/2.0_sp/365.0_sp
                            lfc_ic_h(f_l(1),h_i(i_l,t_l,2),counter_ic_h(f_l(1),h_i(i_l,t_l,2)))=          l_bar(h_i(i_l,t_l,1))/2.0_sp/365.0_sp
                        end if
                        fc_pi_h_it(PI_q_i(i_l),h_i(i_l,t_l,2))=lfc_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2),counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2)))-data_lfc_PI(PI_q_i(i_l),h_i(i_l,t_l,2))
                        fc_ic_h_it(f_l(1),h_i(i_l,t_l,2))=lfc_ic_h(f_l(1),h_i(i_l,t_l,2),counter_ic_h(f_l(1),h_i(i_l,t_l,2)))-data_lfc_IC(f_l(1),h_i(i_l,t_l,2))
                    end if
                    
                    if (h_i(i_l,t_l,2)==1) then
                        counter_ut(1,t_l)=counter_ut(1,t_l)+1
                        assets_ut(1,t_l,counter_ut(1,t_l))=a_it(t_l)
                    else
                        counter_ut(2,t_l)=counter_ut(2,t_l)+1
                        assets_ut(2,t_l,counter_ut(2,t_l))=a_it(t_l)
                    end if
                else 
                    x_it(t_l+1)=-9.0_sp
                    a_it(t_l+1)=-9.0_sp
                end if
                moments_it(:,t_l,i_l)=(/reshape(assets_pi_age_group_it,(/L_PI*obs*groups,1/)),&
                                        reshape(assets_pi_age_group_it_b,(/L_PI*obs*groups,1/)),&
                                        reshape(beq100_it,(/L_PI*f_t,1/)), &
                                        reshape(nw_ic_it,(/f_t*obs*groups,1/)), &
                                        reshape(fc_pi_h_it,(/L_PI*clusters,1/)), &
                                        reshape(fc_ic_h_it,(/f_t*clusters,1/))/)
            end do
        end do
        
        !Compute moments
        do pi_l=1,L_PI; do t_l=1,7; do g_l=1,4
            if (counter_pi_age_group(pi_l,t_l,g_l)>1) then
                call compute_percentile(assets_pi_age_group(pi_l,t_l,g_l,1:counter_pi_age_group(pi_l,t_l,g_l)), &
                                    counter_pi_age_group(pi_l,t_l,g_l),50, & 
                                    med_assets_pi_age_group(pi_l,t_l,g_l,ns))
            end if
            if (counter_pi_age_group_b(pi_l,t_l,g_l)>1) then
                call compute_percentile(assets_pi_age_group_b(pi_l,t_l,g_l,1:counter_pi_age_group_b(pi_l,t_l,g_l)), &
                                    counter_pi_age_group_b(pi_l,t_l,g_l),50, & 
                                    med_assets_pi_age_group_b(pi_l,t_l,g_l,ns))
            end if
        end do; end do; end do
        do pi_l=1,L_PI; do h_l=1,clusters
            if (counter_pi_h(pi_l,h_l)>1) then
                av_lfc_pi_h(pi_l,h_l,ns)=real(sum(lfc_pi_h(pi_l,h_l,1:counter_pi_h(pi_l,h_l))))/real(counter_pi_h(pi_l,h_l))
            end if 
        end do; end do
        do pi_l=1,L_PI;do f_l2=1,f_t
            if (counter_ic(f_l2)>1) then
                av_beq100_ic(pi_l,f_l2,ns)=real(sum(beq100_ic(pi_l,f_l2,1:counter_beq100(pi_l,f_l2))))/real(counter_beq100(pi_l,f_l2))
            end if 
        end do;end do

        do f_l2=1,f_t; do h_l=1,clusters
            if (counter_ic_h(f_l2,h_l)>1) then
                av_lfc_ic_h(f_l2,h_l,ns)=real(sum(lfc_ic_h(f_l2,h_l,1:counter_ic_h(f_l2,h_l))))/real(counter_ic_h(f_l2,h_l))
            end if
        end do; end do
        moments_s(:,ns)=sum(sum(moments_it,3),2)/real(indv)/real(obs)
        do h_l=1,2; do t_l=1,7
            call compute_percentile(assets_ut(h_l,t_l,1:counter_ut(h_l,t_l)), &
                                    counter_ut(h_l,t_l),50, & 
                                    model_NW_h_ut_ns(h_l,t_l,ns))
        end do; end do
        !Wealth moments by ic
        do f_ll=1,f_t;do t_l=1,7;do g_l=1,4
            if (counter_ic_nw(f_ll,t_l,g_l)>1) then
                call compute_percentile(assets_ic(f_ll,t_l,g_l,1:counter_ic_nw(f_ll,t_l,g_l)), &
                                    counter_ic_nw(f_ll,t_l,g_l),75, &
                                    assets_ic_ns(f_ll,t_l,g_l,ns)) 
            end if 
        end do;end do;end do
    end do

    model_moments(:,1)=sum(moments_s,2)/real(samples_per_i)
    open(unit=9,file='model_moments_2.txt')
    do h_l=1,moment_conditions
        write(9,*) model_moments(h_l,1)
    end do
    close(9)
    moments_NW_PI1=sum(med_assets_pi_age_group,4)/real(samples_per_i)
    moments_NW_PI1b=sum(med_assets_pi_age_group_b,4)/real(samples_per_i)
    moments_beq100_IC=sum(av_beq100_ic,3)/real(samples_per_i)
    print*,'moments_beq100_IC, f1',moments_beq100_IC(:,1)
    print*,'moments_beq100_IC, f2',moments_beq100_IC(:,2)
    moments_lfc_PI=sum(av_lfc_pi_h,3)/real(samples_per_i)
    moments_lfc_PI(:,1)=-9.0_sp
    moments_lfc_IC=sum(av_lfc_ic_h,3)/real(samples_per_i)
    moments_lfc_IC(:,1)=-9.0_sp
    moments_NW_IC1=sum(assets_ic_ns,4)/real(samples_per_i)

    model_moments1(:,1)=(/reshape(moments_NW_PI1,(/L_PI*obs*groups,1/)), &
                          reshape(moments_NW_PI1b,(/L_PI*obs*groups,1/)), &
                          reshape(moments_beq100_IC,(/L_PI*f_t,1/)), &
                          reshape(moments_NW_IC1,(/f_t*obs*groups,1/)), &
                          reshape(moments_lfc_PI,(/L_PI*clusters,1/)), &
                          reshape(moments_lfc_IC,(/f_t*clusters,1/))/)
    
    model_NW_h_ut=sum(model_NW_h_ut_ns,3)/real(samples_per_i)
    
    open(unit=9,file='moments_NW_IC1.txt')
        write(9,*) moments_NW_IC1
    close(9)
    
    
end subroutine
    
