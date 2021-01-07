subroutine discretize_var_ep(sigma2_varep,nzz,varep_grid,pr_varep)
    use nrtype
    implicit none
    real(SP),intent(in)::sigma2_varep
    integer,intent(in)::nzz
    real(SP),dimension(nzz,1)::varep_grid,pr_varep
    real(SP),dimension(nzz)::x
    integer::i_l
    real,dimension(nzz-1,1)::dz
    real(SP)::x1,x2,pr_1,pr_2
    
    if (nzz==2) then
        x=(/-0.7071067811,0.7071067811/)
    elseif (nzz==3) then
        x=(/-1.22474487139158,0.0,1.22474487139158/)
    elseif (nzz==4) then
        x=(/-1.650680123,-0.5246476232,0.5246476232,1.650680123/)
    elseif (nzz==5) then
        x=(/-2.02018287,-0.9585724646,0.0,0.9585724646,2.02018287/) 
    elseif (nzz==7)then
        x=(/-2.65196135683523,-1.67355162876747,-0.816287882858964,0.0,0.816287882858964,1.67355162876747,2.65196135683523/)
    elseif(nzz==8)then
        x=(/-2.93063742025724,-1.98165675669584,-1.15719371244678,-0.381186990207322,0.381186990207322,1.15719371244678,1.98165675669584,2.93063742025724/)
    end if
    
    do i_l=1,nzz
        varep_grid(i_l,1)=sqrt(2.0d0)*sqrt(sigma2_varep)*x(i_l)
        if (i_l>1) then
            dz(i_l-1,1)=varep_grid(i_l,1)-varep_grid(i_l-1,1)
        end if
    end do
    
    if (sigma2_varep==0.0_sp)then
        pr_varep=1.0_sp/real(nzz)
    else
        i_l=1
        x1=(varep_grid(i_l,1)+dz(i_l,1)/2.0d0)/sqrt(sigma2_varep)
        pr_varep(i_l,1)=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
    
        i_l=nzz
        x1=(varep_grid(i_l,1)-dz(i_l-1,1)/2.0d0)/sqrt(sigma2_varep)
        pr_varep(i_l,1)=1.0d0-(0.5d0+0.5d0*erf(x1/sqrt(2.0d0)))
    
        do i_l=2,nzz-1
            x1=(varep_grid(i_l,1)-dz(i_l-1,1)/2d0)/sqrt(sigma2_varep)
            pr_1=0.5d0+0.5d0*erf(x1/sqrt(2.0d0))
            x2=(varep_grid(i_l,1)+dz(i_l,1)/2d0)/sqrt(sigma2_varep)
            pr_2=0.5d0+0.5d0*erf(x2/sqrt(2.0d0))
            pr_varep(i_l,1)=pr_2-pr_1
        end do
    end if
    
end subroutine
