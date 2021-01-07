subroutine build_grids
    use grids; use dimensions; use structural_p1;use files_savings
    implicit none
    integer::i_l,j_l
    real(SP)::sigma2_xi
    
    !Grid of cash on hand
    coh_grid(1)=coh_min
    do i_l=2,nkk
        coh_grid(i_l)=(coh_max-coh_min)/real(nkk-1)*real(i_l-1)
    end do
    
    !Upload medical expenditures parameter process
    open(unit=10,file=fdir_med//"med_coef.txt")
        read(10,*) med_coef
    close(10)
    rho=med_coef(13,1)
    sigma2_ep=med_coef(14,1)
    sigma2_ze=med_coef(15,1)
    
    !Grid for persistent and transitory shock of medical expenditures using quadrature methods
    call discretize_shocks(rho,sigma2_ep,sigma2_ze,nzz,ep_grid,pr_ep,xi_grid,pr0_xi,pr_pxi)

end subroutine
        
    