subroutine ECV_V_k_l2(u_x,V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k_l2,V_k2,ECV_k2)
    use dimensions;use nrtype; use structural_p2; use grids; use structural_p1
    implicit none
    integer,intent(in)::ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k_l2
    real(SP),intent(in)::u_x
    real(SP),dimension(nkk,clusters+1,nzz,nzz2,2),intent(in)::V
    real(SP),intent(out)::V_k2,ECV_k2
    integer::h_l2,ps_l2,ts_l2,pos_x2,z_l2
    real(SP)::alpha,x2
    
    ECV_k2=0.0_sp
     do ts_l2=1,nzz;do z_l2=1,nzz2;do ps_l2=1,nzz;do h_l2=1,clusters+1
        if (h_l2<clusters+1) then
            x2=min(max((1.0_sp+r)*coh_grid(k_l2)-m_exp_all(t_l,ge_l,int(real(i_l+1)/2.0_sp),h_l2,ps_l2,ts_l2)+b(i_l,ge_l),0.0_sp),coh_grid(nkk)-0.1_sp)
        else
            x2=coh_grid(k_l2)
        end if
            pos_x2=int(x2/(coh_grid(2)-coh_grid(1))+1.0000001_sp)
            alpha=1.0_sp-(x2-coh_grid(pos_x2))/(coh_grid(2)-coh_grid(1))
            if (coh_grid(pos_x2)==x2) then
                ECV_k2=ECV_k2+pr_ep(ts_l2,1)*pr_pxi(ps_l,ps_l2)*H_av(h_l,h_l2,t_l,int(real(i_l+1)/2.0_sp),ge_l)*V(pos_x2,h_l2,ps_l2,z_l2,2)
            else
                ECV_k2=ECV_k2+pr_varep(z_l2,1)*pr_ep(ts_l2,1)*pr_pxi(ps_l,ps_l2)*H_av(h_l,h_l2,t_l,int(real(i_l+1)/2.0_sp),ge_l)*&
                (alpha*V(pos_x2,h_l2,ps_l2,z_l2,2)+(1.0_sp-alpha)*V(pos_x2+1,h_l2,ps_l2,z_l2,2))
            end if
    end do; end do; end do;end do
    if (x_l>1) then
        V_k2=u_x+beta*ECV_k2
    else
        V_k2=-1.0_sp/0.0_sp
    end if

end subroutine