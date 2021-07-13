subroutine create_tables_draft(parameters)
    use dimensions;use nrtype;use targets;use structural_p1;use structural_p2
    implicit none
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    integer::pi_l,h_l,f_l,m_l,ind
    character(LEN=25),dimension(L_PI)::strings_pi
    character(LEN=18),dimension(clusters)::strings_h
    character(LEN=20),dimension(f_t)::strings_f
    
    real(SP),dimension(f_t,generations-1,2)::p50_75_assets_ic_age
    real(SP),dimension(L_PI,generations-1)::med_assets_pi_age
    real(SP),dimension(L_PI+1)::EDP_bnk,CV_bnk,CV_bnk2,EDP_bnk2
    real(SP),dimension(moment_conditions,1)::model_moments1
    real(SP),dimension(L_PI,clusters)::moments_MD_PI
    real(SP),dimension(f_t,clusters)::moments_lfc_IC,moments_MD_IC
    integer,dimension(3)::v_pi=(/1,3,5/)
    
    strings_pi(1)="Bottom PI quintile"
    strings_pi(2)="Middle PI quintile"
    strings_pi(3)="Top PI quintile"
    strings_h(1)="Healthy"
    strings_h(2)="Physically frail"
    strings_h(3)="Mentally frail "
    strings_h(4)="Impaired"
    strings_f(1)="\textit{Distant}"
    strings_f(2)="\textit{Close}"

    !meeting LTC needs across PI
    open(unit=9,file='C:\Users\jbueren\Google Drive\JMP\Draft\tables\data_LTC_PI.txt')
    write(9,'(A60)'),"Health & \multicolumn{3}{c}{Permanent income quintile} \\  "
    write(9,'(A60)'),"status & Bottom & Middle & Top \\ \hline "
    do h_l=2,clusters
        write(9,'(A18,A2,F5.1,A2,F5.1,A2,F5.1,A4)') strings_h(h_l),"&",data_lfc_PI(v_pi(1),h_l),"&",data_lfc_PI(v_pi(2),h_l),"&",data_lfc_PI(v_pi(3),h_l),"\\"
    end do
    close(9)
    
    !meeting LTC needs across family types
    open(unit=9,file='C:\Users\jbueren\Google Drive\JMP\Draft\tables\data_LTC_IC.txt')
    do f_l=1,f_t;do h_l=2,clusters            
        if (h_l==clusters .and. f_l/=f_t) then
            write(9,'(A20,A2,A18,A2,F5.1,A2,F5.1,A2,F5.1,A2,F5.1,A9)') strings_f(f_l),"&",strings_h(h_l),"&",l_ic(f_l,h_l,1)/365/2,"&",data_lfc_IC(f_l,h_l),"&",pr_nh(f_l,h_l,1,2)*100.0d0,"&",pr_nh(f_l,h_l,2,1)*100.0d0,"\\ [5pt]"
        else
            write(9,'(A20,A2,A18,A2,F5.1,A2,F5.1,A2,F5.1,A2,F5.1,A3)') strings_f(f_l),"&",strings_h(h_l),"&",l_ic(f_l,h_l,1)/365/2,"&",data_lfc_IC(f_l,h_l),"&",pr_nh(f_l,h_l,1,2)*100.0d0,"&",pr_nh(f_l,h_l,2,1)*100.0d0,"\\"
        end if
    end do;end do
    close(9)
    
    !Parameters
    open(unit=9,file='C:\Users\jbueren\Google Drive\JMP\Draft\tables\parameters.txt')
        write(9,'(A50 )') 'Risk Aversion & & \\'
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '$\sigma$: Consumption &', parameters(1),"& (",se(1),")\\"
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '$\nu$: Care hours &', parameters(7),"& (",se(7),")\\ [5pt]"
        write(9,'(A50 )') 'LTC & & \\'
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '$\alpha(h=2)$: Physically frail &', parameters(3),"& (",se(3),")\\"
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '$\alpha(h=3)$: Mentally frail &', parameters(4),"& (",se(4),")\\"
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '$\alpha(h=4)$: Impaired &', parameters(5),"& (",se(5),")\\"
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '$\tau$: Share formal care &', parameters(9),"& (",se(9),")\\"
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '$\omega$: substitution formal/informal care &', parameters(10),"& (",se(10),")\\ [5pt]"
        write(9,'(A50 )') 'Bequest & & \\'
        write(9,'(A50, F9.2, A3, F9.2,A10 )') '$\delta$: curvature &', parameters(6),"& (",se(6),")\\"
        write(9,'(A55, F6.2, A3, F9.2,A10 )') '$\lambda$($F$=\textit{Distant}): marginal utility &', parameters(8),"& (",se(8),")\\ "
        write(9,'(A50, F6.2, A3, F9.2,A10 )') '$\lambda$($F$=\textit{Close}): marginal utility &', parameters(8),"& (",se(8),")\\ [5pt]"
        write(9,'(A70 )') 'Maximum transfer to achieve utility floor $\times 10^3 $ & & \\'
        write(9,'(A50, F6.2, A3, F5.2,A10 )') '\underline{x}\textit{(h=1)}: Healthy &', parameters(2),"& (",se(2),")\\"
    close(9)
    

    
end subroutine