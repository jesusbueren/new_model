subroutine first_step_p
    use structural_p1; use files_savings; use grids
    implicit none
    integer::t_l,ge_l,i_l,h_l,ps_l,ts_l
    
    !Male's mean income across PI q
    b(:,1)=(/6.134_sp,9.504_sp,11.713_sp, 13.826_sp,15.578_sp,17.506_sp, 20.149_sp,23.279_sp,29.254_sp,44.517_sp/)
    !Female's mean income across PI q
    b(:,2)=(/6.055_sp,9.340_sp,11.646_sp, 13.582_sp, 15.543_sp,17.580_sp,20.009_sp,23.360_sp,29.215_sp,42.526_sp/)
    !During two years
    b=b*2.0_sp
      
    !Upload av medical expenditures
    open(unit=10,file=fdir_med//"med_coef.txt")
        read(10,*) med_coef
    close(10)
    
    do t_l=1,generations; do ge_l=1,L_gender;do i_l=1,L_PI2; do h_l=1,clusters; do ps_l=1,nzz; do ts_l=1,nzz
        call med_expenditure(t_l,ge_l,i_l,h_l,xi_grid(ps_l,1),ep_grid(ts_l,1),med_coef,DIM,m_exp_all(t_l,ge_l,i_l,h_l,ps_l,ts_l))
    end do; end do;end do; end do;end do; end do
   
    !Upload health transition probabilities
    open(unit=9,file=path_h//'H_av.txt')
        read(9,*)  H_av
    close(9)
    H_av(:,clusters+1,generations,:,:)=1.0_sp
    H_av(:,1:clusters,generations,:,:)=0.0_sp
    
end subroutine