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
    
   
    
    !Parameters
    open(unit=9,file='C:\Users\jbueren\Google Drive\JMP\Draft\tables\parameters_restricted.txt')
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
    !
    !!Load model moments
    !open(unit=9,file='model_moments.txt')
    !    read(9,*) model_moments1
    !close(9)
    !
    !moments_MD_PI=reshape(model_moments1(L_PI*obs*groups*2+1:L_PI*obs*groups*2+L_PI*clusters,1),(/L_PI,clusters/))
    !moments_MD_IC=reshape(model_moments1(L_PI*obs*groups*2+L_PI*clusters+1:L_PI*obs*groups*2+L_PI*clusters+f_t*clusters,1),(/f_t,clusters/))
    !moments_lfc_IC=reshape(model_moments1(L_PI*obs*groups*2+L_PI*clusters+f_t*clusters+L_PI2*clusters+1:L_PI*obs*groups*2+L_PI*clusters+f_t*clusters+L_PI2*clusters+f_t*clusters,1),(/f_t,clusters/))
    !
    !open(unit=9,file='C:\Users\jbueren\Google Drive\JMP\Draft\tables\model_match.txt')
    !do h_l=1,clusters
    !    do m_l=1,2
    !        if (m_l==1) then
    !            if (h_l==1) then
    !                write(9,'(A13,A3,A3,A3,A3,A3,A3,A2,F5.2,A2,F5.2,A2,F5.2,A2,F5.2,A5)') strings_h(h_l),"&","-","&","-","&","-" &
    !                                                             ,"&",data_MD_PI(1,h_l),"&",data_MD_PI(2,h_l),"&",data_MD_PI(3,h_l),"&",data_MD_PI(4,h_l),"\\"
    !            else
    !                write(9,'(A13,A3,F4.1,A3,F4.1,A3,F4.1,A3,F5.2,A3,F5.2,A3,F5.2,A3,F5.2,A5)') strings_h(h_l),"&",data_lfc_IC(1,h_l),"&",data_lfc_IC(2,h_l),"&",data_lfc_IC(3,h_l) &
    !                                                             ,"&",data_MD_PI(1,h_l),"&",data_MD_PI(2,h_l),"&",data_MD_PI(3,h_l),"&",data_MD_PI(4,h_l),"\\"
    !            end if
    !        else
    !            if (h_l==1) then
    !                write(9,'(A13,A6,A3,A6,A3,A6,A3,A6,F5.2,A6,F5.2,A6,F5.2,A6,F5.2,A10)') "","& [","-"," ] & [","-"," ] & [","-" &
    !                                                             ," ] & [",moments_MD_PI(1,h_l)," ] & [",moments_MD_PI(2,h_l)," ] & [",moments_MD_PI(3,h_l)," ] & [",moments_MD_PI(4,h_l)," ]\\ [5pt]"
    !            else
    !                write(9,'(A13,A6,F4.1,A6,F4.1,A6,F4.1,A6,F5.2,A6,F5.2,A6,F5.2,A6,F5.2,A10)') "","& [",moments_lfc_IC(1,h_l)," ] & [",moments_lfc_IC(2,h_l)," ] & [",moments_lfc_IC(3,h_l) &
    !                                                             ," ] & [",moments_MD_PI(1,h_l)," ] & [",moments_MD_PI(2,h_l)," ] & [",moments_MD_PI(3,h_l)," ] & [",moments_MD_PI(4,h_l)," ]\\ [5pt]"
    !            end if
    !        end if
    !    end do
    !end do
    !close(9)
    
end subroutine