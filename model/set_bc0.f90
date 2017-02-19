module set_bc0_module

   use param1, only: is_defined
   use bl_fort_module, only : c_real
   use iso_c_binding , only: c_int

   contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_bc0                                                 C
!  Purpose: This subroutine does the initial setting of all boundary   C
!  conditions. The user specifications of the boundary conditions are  C
!  checked for veracity in various check_data/ routines:               C
!  (e.g., check_boundary_conditions).                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      subroutine set_bc0(slo, shi, ulo, uhi, vlo, vhi, wlo, whi, &
                         p_g, ep_g, u_g, v_g, w_g, ro_g0, &
                         bc_ilo_type, bc_ihi_type, bc_jlo_type, bc_jhi_type, &
                         bc_klo_type, bc_khi_type, flag)

! Modules
!--------------------------------------------------------------------//
      use bc                , only: bc_u_g, bc_v_g, bc_w_g, bc_p_g, bc_ep_g
      use ic                , only: PINF_, POUT_, MINF_, MOUT_
      use geometry          , only: domlo, domhi

      use scales, only: scale_pressure

      implicit none

      integer(c_int), intent(in   ) :: slo(3),shi(3)
      integer(c_int), intent(in   ) :: ulo(3),uhi(3),vlo(3),vhi(3),wlo(3),whi(3)

      real(c_real), intent(inout) ::  p_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) :: ep_g&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
      real(c_real), intent(inout) ::  u_g&
         (ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
      real(c_real), intent(inout) ::  v_g&
         (vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
      real(c_real), intent(inout) ::  w_g&
         (wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))

      real(c_real), intent(in   ) :: ro_g0

      integer(c_int), intent(in   ) :: bc_ilo_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_ihi_type&
         (slo(2):shi(2),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jlo_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_jhi_type&
         (slo(1):shi(1),slo(3):shi(3),2)
      integer(c_int), intent(in   ) :: bc_klo_type&
         (slo(1):shi(1),slo(2):shi(2),2)
      integer(c_int), intent(in   ) :: bc_khi_type&
         (slo(1):shi(1),slo(2):shi(2),2)
      integer(c_int), intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

! Local variables
!--------------------------------------------------------------------//
! local index for boundary condition
      integer :: bcv, i,j,k

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi

!--------------------------------------------------------------------//

! Incompressible cases require that Ppg specified for one cell.
! The following attempts to pick an appropriate cell.
      CALL SET_IJK_P_G(slo,shi,ro_g0,flag)


      nlft = max(0,domlo(1)-slo(1))
      nbot = max(0,domlo(2)-slo(2))
      ndwn = max(0,domlo(3)-slo(3))

      nrgt = max(0,shi(1)-domhi(1))
      ntop = max(0,shi(2)-domhi(2))
      nup  = max(0,shi(3)-domhi(3))

      if (nlft .gt. 0) then
         ilo = domlo(1)
         do i = 1, nlft
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
                  bcv = bc_ilo_type(j,k,2)
                  if(bc_ilo_type(j,k,1) == PINF_ .or. &
                     bc_ilo_type(j,k,1) == MINF_ .or. &
                     bc_ilo_type(j,k,1) == MOUT_) then

                     p_g(ilo-i,j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(ilo-i,j,k) = bc_ep_g(bcv)

                     u_g(ilo-i,j,k) = bc_u_g(bcv)
                     v_g(ilo-i,j,k) = 0.0d0
                     w_g(ilo-i,j,k) = 0.0d0

                  elseif(bc_ilo_type(j,k,1) == POUT_) then
                     p_g(ilo-i,j,k) = scale_pressure(bc_p_g(bcv))
                  endif
               end do
            end do
         end do
      endif

      if (nrgt .gt. 0) then
         ihi = domhi(1)
         do i = 1, nrgt
            do k=slo(3),shi(3)
               do j=slo(2),shi(2)
                  bcv = bc_ihi_type(j,k,2)
                  if(bc_ihi_type(j,k,1) == PINF_ .or. &
                     bc_ihi_type(j,k,1) == MINF_ .or. &
                     bc_ihi_type(j,k,1) == MOUT_) then

                     p_g(ihi+i,j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(ihi+i,j,k) = bc_ep_g(bcv)

                     u_g(ihi+i-1,j,k) = bc_u_g(bcv)
                     v_g(ihi+i-1,j,k) = 0.0d0
                     w_g(ihi+i-1,j,k) = 0.0d0

                  elseif(bc_ihi_type(j,k,1) == POUT_) then
                     p_g(ihi+i,j,k) = scale_pressure(bc_p_g(bcv))
                  endif
               end do
            end do
         end do
      endif

      if (nbot .gt. 0) then
         jlo = domlo(2)
         do j = 1, nbot
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
                  bcv = bc_jlo_type(i,k,2)
                  if(bc_jlo_type(i,k,1) == PINF_ .or. &
                     bc_jlo_type(i,k,1) == MINF_ .or. &
                     bc_jlo_type(i,k,1) == MOUT_) then

                     p_g(i,jlo-j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,jlo-j,k) = bc_ep_g(bcv)

                     u_g(i,jlo-j,k) = 0.0d0
                     v_g(i,jlo-j,k) = bc_v_g(bcv)
                     w_g(i,jlo-j,k) = 0.0d0

                  elseif(bc_jlo_type(i,k,1) == POUT_) then
                     p_g(i,jlo-j,k) = scale_pressure(bc_p_g(bcv))
                  endif
               end do
            end do
         end do
      endif

      if (ntop .gt. 0) then
         jhi = domhi(2)
         do j = 1, ntop
            do k=slo(3),shi(3)
               do i=slo(1),shi(1)
                  bcv = bc_jhi_type(i,k,2)
                  if(bc_jhi_type(i,k,1) == PINF_ .or. &
                     bc_jhi_type(i,k,1) == MINF_ .or. &
                     bc_jhi_type(i,k,1) == MOUT_) then

                     p_g(i,jhi+j,k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,jhi+j,k) = bc_ep_g(bcv)

                     u_g(i,jhi+j-1,k) = 0.0d0
                     v_g(i,jhi+j-1,k) = bc_v_g(bcv)
                     w_g(i,jhi+j-1,k) = 0.0d0

                  elseif(bc_jhi_type(i,k,1) == POUT_) then
                     p_g(i,jhi+j,k) = scale_pressure(bc_p_g(bcv))
                  endif
               end do
            end do
         end do
      endif

      if (ndwn .gt. 0) then
         klo = domlo(3)
         do k = 1, ndwn
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
                  bcv = bc_klo_type(i,j,2)
                  if(bc_klo_type(i,j,1) == PINF_ .or. &
                     bc_klo_type(i,j,1) == MINF_ .or. &
                     bc_klo_type(i,j,1) == MOUT_) then

                     p_g(i,j,klo-k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,j,klo-k) = bc_ep_g(bcv)

                     u_g(i,j,klo-k) = 0.0d0
                     v_g(i,j,klo-k) = 0.0d0
                     w_g(i,j,klo-k) = bc_w_g(bcv)

                  elseif(bc_klo_type(i,j,1) == POUT_) then
                     p_g(i,j,klo-k) = scale_pressure(bc_p_g(bcv))
                  endif
               end do
            end do
         end do
      endif

      if (nup .gt. 0) then
         khi = domhi(3)
         do k = 1, nup
            do j=slo(2),shi(2)
               do i=slo(1),shi(1)
                  bcv = bc_khi_type(i,j,2)
                  if(bc_khi_type(i,j,1) == PINF_ .or. &
                     bc_khi_type(i,j,1) == MINF_ .or. &
                     bc_khi_type(i,j,1) == MOUT_) then

                     p_g(i,j,khi+k) = scale_pressure(bc_p_g(bcv))
                     ep_g(i,j,khi+k) = bc_ep_g(bcv)

                     u_g(i,j,khi+k-1) = 0.0d0
                     v_g(i,j,khi+k-1) = 0.0d0
                     w_g(i,j,khi+k-1) = bc_w_g(bcv)

                  elseif(bc_khi_type(i,j,1) == POUT_) then
                     p_g(i,j,khi+k) = scale_pressure(bc_p_g(bcv))
                  endif
               end do
            end do
         end do
      endif

   end subroutine set_bc0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_IJK_P_G                                             !
!  Purpose: Pick an appropriate control volume to specify Ppg.         !
!                                                                      !
!  Author: J. Musser                                  Date: 07-Nov-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_IJK_P_G (slo, shi, RO_G0, flag)

      ! IJK location where Ppg is fixed.
      use bc, only: IJK_P_G

      use geometry, only: CYCLIC_X, CYCLIC_X_PD, CYCLIC_X_MF
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD, CYCLIC_Y_MF
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD, CYCLIC_Z_MF
      use geometry, only: domlo, domhi

      use funits, only: DMP_LOG

      use bc, only: BC_DEFINED
      use bc, only: BC_TYPE

      ! MFIX Runtime parameters:
      use param, only: DIMENSION_BC
      use param1, only: UNDEFINED_I, is_undefined

      use compar, only: mype
      use funits, only: UNIT_LOG
      use exit_mod, only: mfix_exit

      implicit none

      integer     , intent(in   ) :: slo(3),shi(3) 
      ! Specified constant gas density.
      real(c_real), intent(in) :: ro_g0
      integer         , intent(in) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)
!--------------------------------------------------------------------//
      INTEGER :: BCV

      CHARACTER(len=7) :: Map
      CHARACTER(len=128) :: lMsg

      INTEGER :: l3
      INTEGER :: l2, u2
      INTEGER :: l1, u1

      LOGICAL, parameter :: setDBG = .FALSE.
      LOGICAL :: dFlag
      INTEGER :: iErr

!--------------------------------------------------------------------//

      dFlag = (DMP_LOG .AND. setDBG)

! Initialize.
      iErr = 0
      IJK_P_G = UNDEFINED_I

! This is not needed for compressible cases.
      IF(IS_UNDEFINED(RO_G0)) THEN
         IF(dFlag) write(*,"(3x,A)")                                   &
            'Compressible: IJK_P_g remaining undefined.'
         return
      ELSEIF(ABS(RO_G0) < EPSILON(0.0d0)) THEN
         IF(dFlag) write(*,"(3x,A)")                                   &
            'No gas phase: IJK_P_g remaining undefined.'
         return
      ENDIF

! If there are no cyclic boundaries, look for a pressure outflow.
      lpBCV: DO BCV = 1, DIMENSION_BC
         IF (.NOT.BC_DEFINED(BCV)) cycle lpBCV
         IF (BC_TYPE(BCV) == 'P_OUTFLOW' .OR. &
             BC_TYPE(BCV) == 'P_INFLOW') THEN
            IF(dFlag) write(*,"(3x,A)")                                &
               'Outflow PC defined: IJK_P_g remaining undefined.'
            RETURN
         ENDIF
      ENDDO lpBCV

! Initialize.
         l3 = UNDEFINED_I
         l2 = UNDEFINED_I; u2=l2
         l1 = UNDEFINED_I; u1=l1

! If there are cyclic boundaries, flag a cell along the positive
! domain extreme in the cyclic direction (e.g., domhi(2)).
      IF(CYCLIC_Y .OR. CYCLIC_Y_PD .OR. CYCLIC_Y_MF) THEN

         Map = 'JKI_MAP'
         l3 = domhi(2)
         l2 = domlo(3);  u2 = domhi(3)
         l1 = domlo(1);  u1 = domhi(1)
         lMsg='Cyclic in Y'

      ELSEIF(CYCLIC_X .OR. CYCLIC_X_PD .OR. CYCLIC_X_MF) THEN

         Map = 'IKJ_MAP'
         l3 = domhi(1)
         l2 = domlo(3);  u2 = domhi(3)
         l1 = domlo(2);  u1 = domhi(2)
         lMsg='Cyclic in X'

      ELSEIF(CYCLIC_Z .OR. CYCLIC_Z_PD .OR. CYCLIC_Z_MF) THEN

         Map = 'KIJ_MAP'
         l3 = domhi(3)
         l2 = domlo(1);  u2 = domhi(1)
         l1 = domlo(2);  u1 = domhi(2)
         lMsg='Cyclic in Z'

      ENDIF

! No cyclic boundaries or pressure outflows. The IJ plane is used in
! this case to maximize search region for 2D problems.
      IF(IS_UNDEFINED(l3)) THEN
         Map = 'KIJ_MAP'
         l3 = max((domhi(3)-domlo(3))/2+1,2)
         l2 = domlo(1);  u2 = domhi(1)
         l1 = domlo(2);  u1 = domhi(2)
         lMsg='Center of domain'
      ENDIF

! Debugging messages.
      IF(dFlag) THEN
         write(*,"(3/,3x,'Map: ',A)") Map
         write(*,"(/5x,'l3:',2x,I4)") l3
         write(*,"( 5x,'l2:',2(2x,I4))") l2, u2
         write(*,"( 5x,'l1:',2(2x,I4))") l1, u1
         write(*,"( 5x,'Msg: ',A)") trim(lMsg)
      ENDIF

! Invoke the search routine.
      CALL IJK_Pg_SEARCH(l3, l2, u2, l1, u1, MAP, dFlag, iErr, flag, slo, shi)

      IF(iErr == 0) RETURN

! Error management.
      IF(DMP_LOG) THEN
         SELECT CASE (iErr)
         CASE ( 1001);  WRITE(UNIT_LOG, 1001); WRITE(*,1001)
         CASE ( 2000);  WRITE(UNIT_LOG, 2000); WRITE(*,2000)
         CASE ( 2001);  WRITE(UNIT_LOG, 2001); WRITE(*,2001)
         CASE ( 2002);  WRITE(UNIT_LOG, 2002); WRITE(*,2002)
         CASE DEFAULT
            WRITE(UNIT_LOG, 1000) iErr
            WRITE(*,1000) iErr
         END SELECT

         WRITE(UNIT_LOG, 9000) MAP(1:1), l3, MAP(2:2),                 &
            l2, u2, MAP(3:3), l1, u1
         WRITE(*, 9000) MAP(1:1), l3, MAP(2:2),                        &
            l2, u2, MAP(3:3), l1, u1

         WRITE(*, 9999)
         WRITE(UNIT_LOG, 9999)

      ENDIF


      CALL MFIX_EXIT(myPE)


 1000 FORMAT(//1X,70('*')/' From: SET_IJK_Pg',/,                       &
         ' Error 1000: Unknown error reported. x', I4.4)

 1001 FORMAT(//1X,70('*')/' From: SET_IJK_Pg',/,                       &
         ' Error 1001: Invalid mapping function.')

 2000 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2000: Unknown error reported from IJK_Pg_SEARCH.')

 2001 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2001: Unable to locate fluid cell in search region.')

 2002 FORMAT(//1X,70('*')/' From: SET_IJK_Pg > IJK_Pg_SEARCH',/,       &
         ' Error 2002: Unable to locate fluid cell owner.')

 9000 FORMAT(/' Search plane information:',/,3x,A1,': ',I8,            &
          2(/3x,A1,': ',I8,' x ',I8))

 9999 FORMAT(/' Fatal Error --> Invoking MFIX_EXIT',/1x,70('*'),2/)

      END SUBROUTINE set_ijk_p_g

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Author: J. Musser                                  Date: 07-Nov-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ijk_pg_search(ll3, ll2, lu2, ll1, lu1, lMAP,          &
         ldFlag, iErr, flag, slo, shi)

! Modules
!--------------------------------------------------------------------//
      ! ijk location where Ppg is fixed.
      use bc, only: IJK_P_g

      use param1, only: undefined_i, is_undefined
      use compar, only: numpes, mype
      implicit none

      integer     , intent(in   ) :: slo(3),shi(3)

      INTEGER, intent(IN   )  :: ll3
      INTEGER, intent(IN   )  :: ll2, lu2
      INTEGER, intent(IN   )  :: ll1, lu1
      LOGICAL, intent(IN   ) :: ldFlag
      INTEGER, intent(  OUT)  :: iErr
      integer, intent(in   ) :: flag&
         (slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),4)

      CHARACTER(len=*), intent(IN) :: lMAP

! Local variables
!--------------------------------------------------------------------//
      INTEGER :: lc2, lS2, lE2
      INTEGER :: lc1, lS1, lE1
      INTEGER :: I, J, K
      LOGICAL :: recheck
      INTEGER :: ijk_pg_owner, proc
      INTEGER :: gijk(0:numPEs-1,3)
      INTEGER :: lpCnt

      CHARACTER(len=32) :: cInt

!--------------------------------------------------------------------//

! Initialize Error Flag
      iErr = 2000

! Initialize the Owner ID
      ijk_pg_owner = UNDEFINED_I

! Set the initial search region, a single cell.
      if(ll1 == lu1) then
         lS1 = ll1
         lE1 = lu1
      else
         lS1 = ll1 + (lu1-ll1)/2 + 1
         lE1 = lS1
      endif

      if(ll2 == lu2) then
         lS2 = ll2
         lE2 = lu2
      else
         lS2 = ll2 + (lu2-ll2)/2 + 1
         lE2 = lS2
      endif

      lpCnt = 1
      recheck = .TRUE.
      do while(recheck)

! Initialize the global ijk array to zero. Resetting this array inside
! this do-loop is most likely overkill. This loop should only cycle
! if gIJK is zero.
         gIJK = undefined_i

! Report debugging information for the search region.
         if(ldFlag) then
            write(*,"(/5x,'Pass: ',I4)") lpCnt
            write(*,"( 5x,'lp2 bounds:',2(2x,I4))")lS2, lE2
            write(*,"( 5x,'lp1 bounds:',2(2x,I4))")lS1, lE1
         endif

         lp2: do lc2 = lS2, lE2
         lp1: do lc1 = lS1, lE1
! Map the loop counters to I/J/K indices.
            SELECT CASE (lMap)
            CASE ('JKI_MAP')
               I=lc1; J=ll3; K=lc2
            CASE ('IKJ_MAP')
               I=ll3; J=lc1; K=lc2
            CASE ('KIJ_MAP')
               I=lc2; J=lc1; K=ll3
            CASE DEFAULT
               iErr = 1001
            END SELECT

! If there is fluid at this location, store the IJK and exit loops.
            if(flag(i,j,k,1)==1) then
               gIJK(myPE,1) = I
               gIJK(myPE,2) = J
               gIJK(myPE,3) = K
               exit lp2
            endif
         enddo lp1
         enddo lp2

! Sync gIJK across all processes. Select the lowest ranked process that
! has gIJK defined. This choice is arbitray and doesn't really matter.
! It just needs to be consistent.
         ! CALL global_all_sum(gIJK)
         proc_lp: do proc=0, numPEs-1
            if(gIJK(proc,1) /= undefined_i) then
               IJK_P_g(1) = gIJK(proc,1)
               IJK_P_g(2) = gIJK(proc,2)
               IJK_P_g(3) = gIJK(proc,3)
               IJK_Pg_Owner = proc
               recheck = .FALSE.
               exit proc_lp
            endif
         enddo proc_lp

! If the proceeding section did not find a fluid cell, expand the search
! area and try again.
         if(recheck) then
            if(lS1 > ll1 .OR. lE1 < lu1 .OR.                           &
               lS2 > ll2 .OR. lE2 < lu2) then
! Expand the 1-axis
               lS1 = max((lS1-1), ll1)
               lE1 = min((lE1+1), lu1)
! Expand the 2-axis
               lS2 = max((lS2-1), ll2)
               lE2 = min((lE2+1), lu2)
! The entire seach plane was checked with no fluid cell identified.
! Force IJK_P_g to undefined for later error checking.
            else
               recheck = .FALSE.
               IJK_P_g = UNDEFINED_I
            endif
         endif
      enddo

! Verify that one fluid cell was detected. Otherwise flag the possible
! errors and return.
      if(IS_UNDEFINED(IJK_P_G(1))) then
         iErr = 2001
         return
      elseif(IS_UNDEFINED(IJK_Pg_Owner)) then
         iErr = 2002
         return
      endif

! If debugging, have PE_IO report some information before the
! data is overwritten.
      if(ldFlag) then
         write(*,"(/3x,'IJK_P_g successfully identified!')")
         cInt=''; write(cInt,*) IJK_Pg_Owner
         write(*,"( 5x,'Owner Rank: ',A)")trim(adjustl(cInt))
         write(*,"(' :: ')", advance='no')
         cInt=''; write(cInt,*) IJK_P_G(1)
         write(*,"('(',A)",advance='no') trim(adjustl(cInt))
         cInt=''; write(cInt,*) IJK_P_G(2)
         write(*,"(',',A)",advance='no') trim(adjustl(cInt))
         cInt=''; write(cInt,*) IJK_P_g(3)
         write(*,"(',',A,')',2/)") trim(adjustl(cInt))
      endif

      IERR = 0

   end subroutine IJK_Pg_SEARCH

   end module set_bc0_module
