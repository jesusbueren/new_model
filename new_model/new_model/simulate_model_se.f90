subroutine simulate_model_se(a_policy,g_policy,lfc_x,beq100_policy, &
                             model_moments,densities)
    use dimensions;use grids; use nrtype; use simulation_input; use structural_p1; use pdfs; use structural_p2;use targets
    implicit none
    integer,dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations),intent(in)::a_policy,g_policy
    real(SP),dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations),intent(in)::beq100_policy
    real(SP),dimension(nkk,clusters,f_t,L_PI2),intent(in)::lfc_x
    real(SP),dimension(moment_conditions,1),intent(out)::densities,model_moments
    !Loop variables
    integer::i_l,t_l,pos_x,xi_l,xi_l2,ts_l2,ind,k2_l,ns,g_l,pi_l,ic_l,h_l,f_l2,nwq_l,s_l,it,f_ll
    integer,dimension(2)::f_l
    !Counter variables
    integer,dimension(L_PI,obs,groups)::counter_pi_age_group,counter_pi_age_group_b
    integer,dimension(L_PI,clusters)::counter_pi_h
    integer,dimension(L_PI,clusters)::counter_pi_h2
    integer,dimension(f_t,clusters)::counter_ic_h
    integer,dimension(f_t,obs,groups)::counter_ic_nw
    !Hours of care variables and medicaid variables
    real(SP),dimension(L_PI,clusters,indv)::lfc_pi_h
    real(SP),dimension(L_PI,clusters,indv)::MD_pi_h
    real(SP),dimension(f_t,clusters,indv)::lfc_ic_h,MD_ic_h
    real(SP),dimension(L_PI,clusters,samples_per_i)::av_MD_pi_h
    real(SP),dimension(L_PI,clusters,samples_per_i)::av_lfc_pi_h
    real(SP),dimension(f_t,clusters,samples_per_i)::av_lfc_ic_h
    real(SP),dimension(f_t,obs,groups,samples_per_i)::assets_ic_ns,pdf_assets_ic_ns
    !Store vector variables for assets
    real(SP),dimension(L_PI,obs,groups,1400)::assets_pi_age_group,assets_pi_age_group_b !1400: maximum number of individuals in a group
    real(SP),dimension(f_t,obs,groups,3000)::assets_ic
    !Model moments
    real(SP),dimension(L_PI,obs,groups)::moments_NW_PI,moments_NW_PI1,moments_NW_PIb,moments_NW_PI1b
    real(SP),dimension(f_t)::moments_NW_IC1
    real(SP),dimension(L_PI,obs,groups,samples_per_i)::med_assets_pi_age_group,fraction_below_pi_age_group,pdf_assets_pi_age_group,&
                                                       med_assets_pi_age_group_b,fraction_below_pi_age_group_b,pdf_assets_pi_age_group_b
    !Moment variables variables
    real(SP),dimension(f_t,obs,groups)::nw_ic_it
    real(SP),dimension(L_PI,clusters)::fc_pi_h_it
    real(SP),dimension(f_t)::beq100_it
    real(SP),dimension(f_t,clusters)::fc_ic_h_it 
    real(SP),dimension(L_PI,obs,groups)::assets_pi_age_group_it,assets_pi_age_group_it_b
    real(SP),dimension(moment_conditions,obs,indv)::moments_it
    real(SP),dimension(moment_conditions,samples_per_i)::moments_s
    !Simulation variables
    real(SP)::u,k2,x_pr_md
    real(SP),dimension(obs+1)::x_it,g_it,a_it
    integer,dimension(1)::seed=254
    character::pause_k
    !Set the seed
    call random_seed(PUT=seed) 
    
    !Initialize value for moments
    moments_s=0.0_sp
    med_assets_pi_age_group=-9.0_sp
    pdf_assets_pi_age_group=-9.0_sp
    med_assets_pi_age_group_b=-9.0_sp
    pdf_assets_pi_age_group_b=-9.0_sp
    av_MD_pi_h=-9.0_sp
    av_lfc_pi_h=-9.0_sp
    av_lfc_ic_h=-9.0_sp
    fraction_below_pi_age_group=0.0_sp
    fraction_below_pi_age_group_b=0.0_sp
    do ns=1,samples_per_i
        !Initialize value for moments in each simulation
        counter_pi_age_group=0
        counter_pi_age_group_b=0
        counter_ic_h=0
        counter_pi_h=0
        counter_pi_h2=0
        counter_ic_nw=0
        assets_pi_age_group=-9.0_sp
        assets_pi_age_group_b=-9.0_sp
        lfc_pi_h=0.0_sp
        lfc_ic_h=0.0_sp
        moments_it=0.0_sp
        MD_pi_h=-9
        MD_iC_h=-9
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
                    f_l(it)=3
                end if
            end do
            !Sample true health and measured health using Kim smoother
            do s_l=1,2;do t_l=obs,1,-1
                if (s_h_i(i_l,t_l,s_l)/=-9.0_sp) then
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
                            if (u<sum(f_h_i(i_l,t_l,1:ind)*H_av(1:ind,h_i(i_l,t_l+1,s_l),generation_i(i_l)+t_l-1,PI_q_i2(i_l),gender_i(i_l)))/sum(H_av(1:clusters,h_i(i_l,t_l+1,1),generation_i(i_l)+t_l-1,PI_q_i2(i_l),gender_i(i_l))*f_h_i(i_l,t_l,:)) .or. ind==clusters) then !Kim smoother
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
            
            g_it=0.0_sp
            !Simulate decisions
            assets_pi_age_group_it=0.0_sp
            assets_pi_age_group_it_b=0.0_sp
            beq100_it=0.0_sp
            nw_ic_it=0.0_sp
            fc_pi_h_it=0.0_sp
            fc_ic_h_it=0.0_sp
            do t_l=1,obs
                beq100_it=0.0_sp
                fc_pi_h_it=0.0_sp
                fc_ic_h_it=0.0_sp
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
                        g_it(t_l)=1
                        k2_l=1
                    else
                        g_it(t_l)=0
                        k2_l=a_policy(pos_x,h_i(i_l,t_l,1),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l(1),generation_i(i_l)+t_l-1)
                    end if
                    if (t_l<obs) then
                        a_it(t_l+1)=coh_grid(k2_l)
                    end if
                    !Medicaid and hours of care
                    counter_ic_h(f_l(1),h_i(i_l,t_l,2))=counter_ic_h(f_l(1),h_i(i_l,t_l,2))+1
                    counter_pi_h2(PI_q_i(i_l),h_i(i_l,t_l,2))=counter_pi_h2(PI_q_i(i_l),h_i(i_l,t_l,2))+1
                    counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2))=counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2))+1
                
                    if (g_it(t_l)==0) then
                        lfc_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2),counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2)))=lfc_x(pos_x-k2_l+1,h_i(i_l,t_l,1),f_l(1),PI_q_i2(i_l))/2.0_sp/365.0_sp
                        lfc_ic_h(f_l(1),h_i(i_l,t_l,2),counter_ic_h(f_l(1),h_i(i_l,t_l,2)))=          lfc_x(pos_x-k2_l+1,h_i(i_l,t_l,1),f_l(1),PI_q_i2(i_l))/2.0_sp/365.0_sp
                    else
                        lfc_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2),counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2)))=l_bar(h_i(i_l,t_l,1))/2.0_sp/365.0_sp
                        lfc_ic_h(f_l(1),h_i(i_l,t_l,2),counter_ic_h(f_l(1),h_i(i_l,t_l,2)))=          l_bar(h_i(i_l,t_l,1))/2.0_sp/365.0_sp
                    end if
                    !Moments
                    !Store Assets for computing moments
                    !if (iwendy_i(i_l)==1998 .or. iwendy_i(i_l)==1999) then
                        counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l))=counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l))+1
                        assets_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l),counter_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l)))=a_it(t_l)    
                        if (a_it(t_l)<=data_NW_PI1(PI_q_i(i_l),t_l,group_i(i_l))) then
                            fraction_below_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l),ns)=fraction_below_pi_age_group(PI_q_i(i_l),t_l,group_i(i_l),ns)+1.0_sp
                            assets_pi_age_group_it(PI_q_i(i_l),t_l,group_i(i_l))=0.5_sp
                        else
                            assets_pi_age_group_it(PI_q_i(i_l),t_l,group_i(i_l))=-0.5_sp
                        end if
                    !else
                    !    counter_pi_age_group_b(PI_q_i(i_l),t_l,group_i(i_l))=counter_pi_age_group_b(PI_q_i(i_l),t_l,group_i(i_l))+1
                    !    assets_pi_age_group_b(PI_q_i(i_l),t_l,group_i(i_l),counter_pi_age_group_b(PI_q_i(i_l),t_l,group_i(i_l)))=a_it(t_l)    
                    !    if (a_it(t_l)<=data_NW_PI1b(PI_q_i(i_l),t_l,group_i(i_l))) then
                    !        fraction_below_pi_age_group_b(PI_q_i(i_l),t_l,group_i(i_l),ns)=fraction_below_pi_age_group_b(PI_q_i(i_l),t_l,group_i(i_l),ns)+1.0_sp
                    !        assets_pi_age_group_it_b(PI_q_i(i_l),t_l,group_i(i_l))=0.5_sp
                    !    else
                    !        assets_pi_age_group_it_b(PI_q_i(i_l),t_l,group_i(i_l))=-0.5_sp
                    !    end if
                    !end if
   
                    counter_ic_nw(f_l(2),t_l,group_i(i_l))=counter_ic_nw(f_l(2),t_l,group_i(i_l))+1
                    assets_ic(f_l(2),t_l,group_i(i_l),counter_ic_nw(f_l(2),t_l,group_i(i_l)))=a_it(t_l)
                    if (a_it(t_l)<=data_NW_IC1(f_l(2),t_l,group_i(i_l))) then
                        nw_ic_it(f_l(2),t_l,group_i(i_l))=0.25_sp
                    else
                        nw_ic_it(f_l(2),t_l,group_i(i_l))=-0.75_sp
                    end if 
 
                    fc_pi_h_it(PI_q_i(i_l),h_i(i_l,t_l,2))=lfc_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2),counter_pi_h(PI_q_i(i_l),h_i(i_l,t_l,2)))-data_lfc_PI(PI_q_i(i_l),h_i(i_l,t_l,2))
                    fc_ic_h_it(f_l(1),h_i(i_l,t_l,2))=lfc_ic_h(f_l(1),h_i(i_l,t_l,2),counter_ic_h(f_l(1),h_i(i_l,t_l,2)))-data_lfc_IC(f_l(1),h_i(i_l,t_l,2))
                    beq100_it(f_l(2))=beq100_policy(pos_x,h_i(i_l,t_l,1),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l(1),generation_i(i_l)+t_l-1)-data_beq100_IC(f_l(2))
                else 
                    x_it(t_l+1)=-9.0_sp
                    a_it(t_l+1)=-9.0_sp
                end if
                moments_it(:,t_l,i_l)=(/reshape(assets_pi_age_group_it,(/L_PI*obs*groups,1/)),&
                                        reshape(assets_pi_age_group_it_b,(/L_PI*obs*groups,1/)),&
                                        reshape(beq100_it,(/f_t,1/)), &
                                        reshape(nw_ic_it,(/f_t*obs*groups,1/)), &
                                        reshape(fc_pi_h_it,(/L_PI*clusters,1/)), &
                                        reshape(fc_ic_h_it,(/f_t*clusters,1/))/)
            end do
        end do
        
        !Compute moments
        fraction_below_pi_age_group(:,:,:,ns)=fraction_below_pi_age_group(:,:,:,ns)/real(counter_pi_age_group)
        do pi_l=1,L_PI; do t_l=1,7; do g_l=1,4
            if (counter_pi_age_group(pi_l,t_l,g_l)>1) then
                call compute_percentile(assets_pi_age_group(pi_l,t_l,g_l,1:counter_pi_age_group(pi_l,t_l,g_l)), &
                                        counter_pi_age_group(pi_l,t_l,g_l),50, & 
                                        med_assets_pi_age_group(pi_l,t_l,g_l,ns))
                if (data_NW_PI1(pi_l,t_l,g_l)/=-9.0_sp)then
                    call compute_pdf(assets_pi_age_group(pi_l,t_l,g_l,1:counter_pi_age_group(pi_l,t_l,g_l)),counter_pi_age_group(pi_l,t_l,g_l),med_assets_pi_age_group(pi_l,t_l,g_l,ns),pdf_assets_pi_age_group(pi_l,t_l,g_l,ns))
                end if
                med_assets_pi_age_group(pi_l,t_l,g_l,ns)=med_assets_pi_age_group(pi_l,t_l,g_l,ns)*counter_pi_age_group(pi_l,t_l,g_l)/real(obs)/real(indv)
            end if
            !if (counter_pi_age_group_b(pi_l,t_l,g_l)>1) then
            !    call compute_percentile(assets_pi_age_group_b(pi_l,t_l,g_l,1:counter_pi_age_group_b(pi_l,t_l,g_l)), &
            !                            counter_pi_age_group_b(pi_l,t_l,g_l),50, & 
            !                            med_assets_pi_age_group_b(pi_l,t_l,g_l,ns))
            !    if (data_NW_PI1b(pi_l,t_l,g_l)/=-9.0_sp)then
            !        call compute_pdf(assets_pi_age_group_b(pi_l,t_l,g_l,1:counter_pi_age_group_b(pi_l,t_l,g_l)),counter_pi_age_group_b(pi_l,t_l,g_l),med_assets_pi_age_group_b(pi_l,t_l,g_l,ns),pdf_assets_pi_age_group_b(pi_l,t_l,g_l,ns))
            !    end if
            !    med_assets_pi_age_group_b(pi_l,t_l,g_l,ns)=med_assets_pi_age_group_b(pi_l,t_l,g_l,ns)*counter_pi_age_group_b(pi_l,t_l,g_l)/real(obs)/real(indv)
            !end if
        end do; end do; end do
        !do f_ll=1,f_t;do t_l=1,7;do g_l=1,4
        !    call compute_percentile(assets_ic(f_ll,t_l,g_l,1:counter_ic_nw(f_ll,t_l,g_l)), &
        !                                counter_ic_nw(f_ll,t_l,g_l),75, & 
        !                                assets_ic_ns(f_ll,t_l,g_l,ns))
        !    call compute_pdf(assets_ic(f_ll,t_l,g_l,1:counter_ic_nw(f_ll,t_l,g_l)),counter_ic_nw(f_ll,t_l,g_l),assets_ic_ns(f_ll,t_l,g_l,ns),pdf_assets_ic_ns(f_ll,t_l,g_l,ns))
        !end do;end do;end do
        do pi_l=1,L_PI; do h_l=1,clusters
            if (counter_pi_h(pi_l,h_l)>1) then
                av_lfc_pi_h(pi_l,h_l,ns)=real(sum(lfc_pi_h(pi_l,h_l,1:counter_pi_h(pi_l,h_l))))/real(counter_pi_h(pi_l,h_l))
            end if 
        end do; end do
        do pi_l=1,L_PI; do h_l=1,clusters
            if (counter_pi_h2(pi_l,h_l)>1) then
                av_MD_pi_h(pi_l,h_l,ns)=real(sum(MD_pi_h(pi_l,h_l,1:counter_pi_h(pi_l,h_l))))/real(counter_pi_h(pi_l,h_l))
            end if 
        end do; end do
        
        do f_l2=1,f_t; do h_l=1,clusters
            if (counter_ic_h(f_l2,h_l)>1) then
                av_lfc_ic_h(f_l2,h_l,ns)=real(sum(lfc_ic_h(f_l2,h_l,1:counter_ic_h(f_l2,h_l))))/real(counter_ic_h(f_l2,h_l))
            end if
        end do; end do
        moments_s(:,ns)=sum(sum(moments_it,3),2)/real(indv)/real(obs)
    end do

    model_moments(:,1)=sum(moments_s,2)/real(samples_per_i)
    moments_NW_PI1=sum(med_assets_pi_age_group,4)/real(samples_per_i)
    moments_NW_PI1b=sum(med_assets_pi_age_group_b,4)/real(samples_per_i)

    model_moments(1:2*L_PI*obs*groups,1)=(/ reshape(moments_NW_PI1,(/L_PI*obs*groups,1/)),&
                                            reshape(moments_NW_PI1b,(/L_PI*obs*groups,1/)) /)
    densities=1.0_sp
    densities(1:L_PI*obs*groups*2,1)=(/reshape(sum(pdf_assets_pi_age_group,4)/real(samples_per_i),(/L_PI*obs*groups,1/)),&
                                       reshape(sum(pdf_assets_pi_age_group_b,4)/real(samples_per_i),(/L_PI*obs*groups,1/)) /)  
    
    densities(L_PI*obs*groups*2+f_t+1:L_PI*obs*groups*2+f_t+f_t*obs*groups,1)=(/reshape(sum(pdf_assets_ic_ns,4)/real(samples_per_i),(/f_t*obs*groups,1/))/)
    
    
end subroutine
    
