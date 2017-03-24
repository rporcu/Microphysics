!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DISPLAY_RESID(NIT, IER)                                !
!  Purpose: Display residuals                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DISPLAY_RESID(NIT, resid)&
         bind(C, name="display_resid")

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding, only: c_int

      use residual, only: group_resid

      implicit none

! iteration number
      integer(c_int), intent(in) :: nit
      real(c_real),   intent(in) :: resid(8,2)

      if(group_resid) then
         call display_group_resid
      else
         call display_field_resid
      endif

      return

   contains


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DISPLAY_FIELD_RESID(NIT, IER)                          !
!  Purpose: Display residuals for each field variable.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DISPLAY_FIELD_RESID

      use param1, only: undefined_i
      use residual, only: resid_string, resid_index

      use error_manager, only: err_msg, flush_err_msg

      implicit none

      integer :: ll, lc, ls, le

      if(nit == 1) then
         write(err_msg(1)(1:5),'("  Nit")')
         lc=1
         do ll = 1, 8
            if (resid_index(ll,1) /= undefined_i) then
               ls= 6+10*(lc-1)
               le= 5+10*(lc)
               write(err_msg(1)(ls:le),'(5x,a4)') resid_string(ll)
               lc=lc+1
            endif
         end do
         if(resid_index(8,1) == undefined_i) then
            ls= 6+10*(lc-1)
            le= 5+10*(lc)
            write(err_msg(1)(ls:le),'(2x,a7)') 'Max res'
         endif
         call flush_err_msg(header=.false., footer=.false., log=.false.)
      endif


      write(err_msg(1)(1:5),'(i5)') nit
      lc=1
      do ll = 1, 8
         if(resid_index(ll,1) /= undefined_i) then
            ls= 6+10*(lc-1)
            le= 5+10*(lc)
            write(err_msg(1)(ls:le),'(2x,1pg8.1)') &
               resid(resid_index(ll,1),1)
            lc=lc+1
         endif
      enddo
      if(resid_index(8,1) == undefined_i) then
         ls= 6+10*(lc-1)
         le= 3+10*(lc)
         write(err_msg(1)(ls:le),'(4x,a4)') resid_string(8)
      endif
      call flush_err_msg(header=.false., footer=.false., log=.false.)

      RETURN

      END SUBROUTINE DISPLAY_FIELD_RESID



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DISPLAY_GROUP_RESID(NIT, IER)                          !
!  Author: M. Syamlal                                 Date: 8-JUL-96   !
!                                                                      !
!  Purpose: Display residuals grouped by equation type.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DISPLAY_GROUP_RESID

      use residual, only: RESID_STRING
      use residual, only: RESID_GRP, RESID_GRP_STRING
      use residual, only: HYDRO_GRP

      use error_manager, only: finl_err_msg, err_msg, flush_err_msg, init_err_msg, ivar

      IMPLICIT NONE


      INTEGER :: LC, LS, LE

      if (nit == 1) then
         write(err_msg(1)(1:5),'("  Nit")')
         lc=1

         ls= 6+10*(lc-1)
         le= 5+10*(lc)
         write(err_msg(1)(ls:le),1000) resid_grp_string(hydro_grp)
         lc=lc+1

         ls= 6+10*(lc-1)
         le= 5+10*(lc)
         write(err_msg(1)(ls:le),1000) 'Max res '

         call flush_err_msg(header=.false., footer=.false., log=.false.)
      endif

 1000 format(3x,a7)



      write(err_msg(1)(1:5),'(i5)') nit
      lc=1

      ls= 6+10*(lc-1)
      le= 5+10*(lc)
      write(err_msg(1)(ls:le),1100) resid_grp(hydro_grp)
      lc=lc+1


      ls= 6+10*(lc-1)
      le= 3+10*(lc)
      write(err_msg(1)(ls:le),'(4x,a4)') resid_string(8)

 1100 format(2x,1pg8.1)

      call flush_err_msg(header=.false., footer=.false., log=.false.)

      RETURN
      END SUBROUTINE DISPLAY_GROUP_RESID

      END SUBROUTINE DISPLAY_RESID
