!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: desgrid                                                     !
!  Author: Pradeep G                                                   !
!                                                                      !
!  Purpose: Defines the desgrid and sets the indices; also sets the    !
!  communication between the desgrid for ghost cell exchange.          !
!                                                                      !
!  Comments: The variables are named similar to the fluid grid. More   !
!  descriptions about naming can be found in compar_mod under          !
!  dmp_modules folder.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DESGRID

      use compar, only: iend1, jend1, kend1
      use compar, only: iend1_all, jend1_all, kend1_all
      use compar, only: iend2, jend2, kend2
      use compar, only: istart1_all, jstart1_all, kstart1_all
      use compar, only: istart2, jstart2, kstart2
      use compar, only: nodesi, nodesj, nodesk, numpes, mype
      use des_allocate, only: add_pair
      use discretelement, only: DESGRIDSEARCH_IMAX
      use discretelement, only: DESGRIDSEARCH_JMAX
      use discretelement, only: DESGRIDSEARCH_KMAX
      use discretelement, only: DES_PERIODIC_WALLS_X
      use discretelement, only: DES_PERIODIC_WALLS_Y
      use discretelement, only: DES_PERIODIC_WALLS_Z
      use discretelement, only: DIMN, factor_RLM, max_pip, max_isize, pip
      use discretelement, only: XE, YN, ZT
      use discretelement, only: entering_ghost, normal_ghost, normal_particle
      use discretelement, only: nonexistent, exiting_ghost, entering_particle
      use error_manager, only: err_msg, flush_err_msg, init_err_msg, finl_err_msg
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH
      use geometry, only: imax2, jmax2, kmax2
      use geometry, only: imin1, jmin1, kmin1
      use geometry, only: imin2, kmin2

      implicit none



      contains


!------------------------------------------------------------------------
! subroutine       : desgrid_neigh_build ()
! Purpose          : This particles build the neigh list for the particles
!                    currently active in the system
!------------------------------------------------------------------------
      subroutine desgrid_neigh_build(des_pos_new, pijk, pinc, particle_state,&
         des_radius, neighbor_index)

         use geometry, only: dx, dy, dz
         use discretelement, only: pic
      implicit none

      double precision, intent(in) :: des_radius(:)
      double precision, intent(in) :: des_pos_new(:,:)
      integer, intent(inout) :: neighbor_index(:)
      integer, intent(in) :: particle_state(:)
      integer, intent(in) :: pijk(:,:)
      integer, intent(in) :: pinc(:,:,:)

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lcurpar
      integer lijk,lic,ljc,lkc,li,lj,lk,ltotpic,lpicloc,lneigh,cc
      double precision lsearch_rad,ldistsquared
      double precision :: ldistvec(3)
      double precision :: lcurpar_pos(3)
      double precision :: lcur_off
      integer il_off,iu_off,jl_off,ju_off,kl_off,ku_off
      integer, dimension(:), allocatable :: tmp_neigh
      double precision :: Oodx, Oody, Oodz
!-----------------------------------------------

! loop through neighbours and build the contact particles list for particles
! present in the system

      allocate(tmp_neigh(max_pip*max_pip))

      Oodx = 1.0d0/dx
      Oody = 1.0d0/dy
      Oodz = 1.0d0/dz

      do lcurpar =1,max_pip

         if (neighbor_index(lcurpar) .eq. 0) then
            if (lcurpar .eq. 1) then
               neighbor_index(lcurpar) = 1
            else
               neighbor_index(lcurpar) = neighbor_index(lcurpar-1)
            endif
         endif

         if (nonexistent==particle_state(lcurpar) .or.&
            entering_particle==particle_state(lcurpar) .or. &
            entering_ghost==particle_state(lcurpar) .or.  &
            normal_ghost==particle_state(lcurpar) .or. &
            exiting_ghost==particle_state(lcurpar)) cycle

         lic = PIJK(lcurpar,1)
         ljc = PIJK(lcurpar,2)
         lkc = PIJK(lcurpar,3)

         il_off = 1
         iu_off = 1
         jl_off = 1
         ju_off = 1
         kl_off = 1
         ku_off = 1

         lcurpar_pos(:) = des_pos_new(lcurpar,:)
!   The desgrid size should not be less than 2*dia*rlm_factor
         lcur_off = (lcurpar_pos(1))*Oodx - &
            floor((lcurpar_pos(1))*Oodx)
         if(lcur_off .ge. 0.5) then
            il_off = 0
         else
            iu_off = 0
         endif

         lcur_off = (lcurpar_pos(2))*Oody - &
            floor((lcurpar_pos(2))*Oody)
         if(lcur_off .ge. 0.5) then
            jl_off = 0
         else
            ju_off = 0
         endif

         lcur_off = (lcurpar_pos(3))*Oodz - &
            floor((lcurpar_pos(3))*Oodz)
         if(lcur_off .ge. 0.5) then
            kl_off = 0
         else
            ku_off = 0
         endif

         do lk = lkc-kl_off,lkc+ku_off
         do lj = ljc-jl_off,ljc+ju_off
         do li = lic-il_off,lic+iu_off
            ltotpic =pinc(li,lj,lk)
            do lpicloc = 1,ltotpic
               lneigh = pic(li,lj,lk)%p(lpicloc)
               tmp_neigh(lpicloc) = 0
! Only skip real particles otherwise collisions with ghost, entering,
! and exiting particles are missed.
               if (lneigh .eq. lcurpar) cycle
               if (lneigh < lcurpar .and.normal_particle==particle_state(lneigh)) cycle
               if (nonexistent==particle_state(lneigh)) THEN
                  cycle
               endif

               lsearch_rad = factor_RLM*(des_radius(lcurpar)+des_radius(lneigh))
               ldistvec(1) = lcurpar_pos(1)-des_pos_new(lneigh,1)
               ldistvec(2) = lcurpar_pos(2)-des_pos_new(lneigh,2)
               ldistvec(3) = lcurpar_pos(3)-des_pos_new(lneigh,3)
               ldistsquared = dot_product(ldistvec,ldistvec)
               if (ldistsquared.gt.lsearch_rad*lsearch_rad) cycle
               tmp_neigh(lpicloc) = lneigh
            end do

            do lpicloc = 1,ltotpic
               lneigh = tmp_neigh(lpicloc)
               if (lneigh .ne. 0) then
                  cc = add_pair(lcurpar, lneigh)
               endif
            end do
         end do
         end do
         end do
      end do

      deallocate(tmp_neigh)

   end subroutine desgrid_neigh_build

end module DESGRID
