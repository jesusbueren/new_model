subroutine std_vfi(u_x,k2_l_min,x_l,ge_l,i_l,f_l,t_l,h_l,ps_l,V,V_wo_MD,k2_wo_MD,beq_aux,beq_wo_md)
    use dimensions; use nrtype
    implicit none
    integer,intent(in)::k2_l_min,x_l,ge_l,i_l,f_l,t_l,h_l,ps_l
    real(SP),dimension(nkk,clusters+1,nzz,2),intent(in)::V,beq_aux
    real(SP),dimension(nkk,clusters,f_t,L_PI2),intent(in)::u_x
    real(SP),intent(out)::V_wo_MD,beq_wo_md
    integer,intent(out)::k2_wo_MD
    real(SP),dimension(nkk)::ECV,V_k2,u,beq_k2
    integer::k2_l_max,k2_l
    
        k2_l_max=x_l
        u=-9.0_sp
        ECV=-9.0_sp
        V_k2=-9.0_sp        
        do k2_l=k2_l_min,k2_l_max
            u(k2_l)=u_x(x_l-k2_l+1,h_l,f_l,i_l)
            call ECV_V_k_l2(u(k2_l),V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k2_l,V_k2(k2_l),ECV(k2_l),beq_aux,beq_k2(k2_l))
        end do 
        V_wo_MD=maxval(V_k2(k2_l_min:k2_l_max))
        k2_wo_MD=maxloc(V_k2(k2_l_min:k2_l_max),1)+k2_l_min-1
        beq_wo_md=beq_k2(k2_wo_MD)
        

end subroutine    