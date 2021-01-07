subroutine std_vfi(u_x,k2_l_min,x_l,ge_l,i_l,f_l,t_l,h_l,ps_l,z_l,V,V_wo_MD,k2_wo_MD)
    use dimensions; use nrtype
    implicit none
    integer,intent(in)::k2_l_min,x_l,ge_l,i_l,f_l,t_l,h_l,ps_l,z_l
    real(SP),dimension(nkk,clusters+1,nzz,nzz2,2),intent(in)::V
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2),intent(in)::u_x
    real(SP),intent(out)::V_wo_MD
    integer,intent(out)::k2_wo_MD
    real(SP),dimension(nkk)::ECV,V_k2,u
    integer::k2_l_max,k2_l
    
        k2_l_max=x_l
        u=-9.0_sp
        ECV=-9.0_sp
        V_k2=-9.0_sp        
        do k2_l=k2_l_min,k2_l_max
            u(k2_l)=u_x(x_l-k2_l+1,h_l,z_l,f_l,i_l)
            call ECV_V_k_l2(u(k2_l),V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k2_l,V_k2(k2_l),ECV(k2_l))
        end do
        V_wo_MD=maxval(V_k2(k2_l_min:k2_l_max))
        k2_wo_MD=maxloc(V_k2(k2_l_min:k2_l_max),1)+k2_l_min-1

end subroutine    