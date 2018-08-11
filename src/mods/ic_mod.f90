!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ic                                                          !
!  Purpose: Global initial conditions variables.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
module ic

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding , only: c_int

  use param, only: dim_ic, dim_m, dim_n_g, dim_n_s


  ! Boundary condition coordinates
  real(rt) :: IC_X_w(dim_ic), IC_X_e(dim_ic)
  real(rt) :: IC_Y_s(dim_ic), IC_Y_n(dim_ic)
  real(rt) :: IC_Z_b(dim_ic), IC_Z_t(dim_ic)

  ! Void fraction in a specified boundary
  real(rt) :: IC_EP_g(dim_ic), IC_EP_s(dim_ic, dim_m)

  ! Initial gas pressure
  real(rt) :: IC_P_g(dim_ic)

  ! Initial velocities in specified region
  real(rt) :: IC_U_g(dim_ic), IC_U_s(dim_ic, dim_m)
  real(rt) :: IC_V_g(dim_ic), IC_V_s(dim_ic, dim_m)
  real(rt) :: IC_W_g(dim_ic), IC_W_s(dim_ic, dim_m)

  ! Heat transfer boundary condition
  real(rt) :: IC_T_g(dim_ic), IC_T_s(dim_ic, dim_m)

  ! Species transfer boundary condition
  real(rt) :: IC_X_g(dim_ic, dim_n_g), IC_X_s(dim_ic, dim_m, dim_n_s)

  ! Particle Size properties
  character(len=16) :: ic_dp_dist(dim_ic, dim_m)
  real(rt) :: ic_dp_mean(dim_ic, dim_m)
  real(rt) :: ic_dp_std(dim_ic, dim_m)
  real(rt) :: ic_dp_min(dim_ic, dim_m)
  real(rt) :: ic_dp_max(dim_ic, dim_m)

  ! Particle density properties
  character(len=16) :: ic_ro_s_dist(dim_ic, dim_m)
  real(rt) :: ic_ro_s_mean(dim_ic, dim_m)
  real(rt) :: ic_ro_s_std(dim_ic, dim_m)
  real(rt) :: ic_ro_s_min(dim_ic, dim_m)
  real(rt) :: ic_ro_s_max(dim_ic, dim_m)

  character(len=16) :: ic_pack(dim_ic)

contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: ic_defined                                               !
!                                                                      !
! Purpose: Return if a IC region has been defined based on coordinates !
! defined in the input deck.                                           !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
  logical function ic_defined(icv)

    use param, only: is_defined

    integer, intent(in) :: icv

    ic_defined = &
         is_defined(ic_x_w(icv)) .or. is_defined(ic_x_e(icv)) .or. &
         is_defined(ic_y_s(icv)) .or. is_defined(ic_y_n(icv)) .or. &
         is_defined(ic_z_b(icv)) .or. is_defined(ic_z_t(icv))

   end function ic_defined

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: write_out_ic                                            !
!                                                                      !
!  Purpose: Echo user input for IC regions.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   subroutine write_out_ic(unit_out, dx, dy, dz)

     use calc_cell_module, only: calc_cell_ic
     use param, only: zero, is_defined
     use fld_const, only: ro_g0

     integer,        intent(in) :: unit_out
     real(rt)  , intent(in) :: dx, dy, dz

     integer :: icv, m
     integer :: i_w, j_s, k_b
     integer :: i_e, j_n, k_t

     write (unit_out, 1500)
1500 format(//,3X,'6. INITIAL CONDITIONS')

     do icv = 1, dim_ic
        if (ic_defined(icv)) then

           write (unit_out, 1510) icv
1510  format(/7X,'Initial condition no : ',I4)

           call calc_cell_ic(dx, dy, dz, &
             ic_x_w(icv), ic_y_s(icv), ic_z_b(icv), &
             ic_x_e(icv), ic_y_n(icv), ic_z_t(icv), &
             i_w, i_e, j_s, j_n, k_b, k_t)

            write (unit_out, 1520) &
               ic_x_w(icv), dx*dble(i_w), ic_x_e(icv), dx*dble(i_e+1), &
               ic_y_s(icv), dy*dble(j_s), ic_y_n(icv), dy*dble(j_n+1), &
               ic_z_b(icv), dz*dble(k_b), ic_z_t(icv), dz*dble(k_t+1)

1520 format(9x,45X,' Specified  ',5X,' Simulated  ',/&
         9X,'X coordinate of west face   (IC_X_w) ...... ',g12.5, 5x, g12.5/,&
         9x,'X coordinate of east face   (IC_X_e) ...... ',g12.5, 5x, g12.5/,&
         9x,'Y coordinate of south face  (IC_Y_s) ...... ',g12.5, 5x, g12.5/,&
         9x,'Y coordinate of north face  (IC_Y_n) ...... ',g12.5, 5x, g12.5/,&
         9x,'Z coordinate of bottom face (IC_Z_b) ...... ',g12.5, 5x, g12.5/,&
         9x,'Z coordinate of top face    (IC_Z_t) ...... ',g12.5, 5x, g12.5/)

            write (unit_out, 1530) i_w, i_e, j_s, j_n, k_b, k_t

1530  format(&
         9X,'I index of cell at west   (IC_I_w) ',24('.'),1x,I4,/,&
         9X,'I index of cell at east   (IC_I_e) ',24('.'),1x,I4,/,&
         9X,'J index of cell at south  (IC_J_s) ',24('.'),1x,I4,/,&
         9X,'J index of cell at north  (IC_J_n) ',24('.'),1x,I4,/,&
         9X,'K index of cell at bottom (IC_K_b) ',24('.'),1x,I4,/,&
         9X,'K index of cell at top    (IC_K_t) ',24('.'),1x,I4)

            if(ro_g0 > zero) then
               write(unit_out, "(' ')")
               write (unit_out, 1640) ic_ep_g(icv)
               write(unit_out, "(' ')")
               if(is_defined(ic_p_g(icv))) &
                 write (unit_out, 1641) ic_p_g(icv)
               write (unit_out, 1650) ic_u_g(icv)
               write (unit_out, 1651) ic_v_g(icv)
               write (unit_out, 1652) ic_w_g(icv)
            endif
1640  format(9X,'Gas phase volume fraction (IC_EP_g) ....... ',g12.5)
1641  format(9X,'Gas pressure (IC_P_g) ..................... ',g12.5)
1650  format(9X,'X-component of gas velocity (IC_U_g) ...... ',g12.5)
1651  format(9X,'Y-component of gas velocity (IC_V_g) ...... ',g12.5)
1652  format(9X,'Z-component of gas velocity (IC_W_g) ...... ',g12.5)

            do m = 1, dim_m
               if(ic_ep_s(icv,m) > zero) then
                  write(unit_out, "(' ')")
                  write(unit_out, 1660) m, ic_ep_s(icv,m)
                  write(unit_out, "(' ')")
                  write(unit_out,1670)m,ic_u_s(icv,m)
                  write(unit_out,1671)m,ic_v_s(icv,m)
                  write(unit_out,1672)m,ic_w_s(icv,m)
               endif
            enddo

1660  format(9X,'Solids phase-',I2,' Volume fraction (IC_EP_s) ............. ',g12.5)
1670  format(9X,'X-component of solids phase-',I2,' velocity (IC_U_s) ...... ',g12.5)
1671  format(9X,'Y-component of solids phase-',I2,' velocity (IC_V_s) ...... ',g12.5)
1672  format(9X,'Z-component of solids phase-',I2,' velocity (IC_W_s) ...... ',g12.5)

         endif
      enddo

    end subroutine write_out_ic

end module ic
