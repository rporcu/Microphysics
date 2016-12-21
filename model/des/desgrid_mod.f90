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
      use discretelement, only: dg_pic
      use discretelement, only: entering_ghost, normal_ghost, normal_particle
      use discretelement, only: nonexistent, exiting_ghost, entering_particle
      use error_manager, only: err_msg, flush_err_msg, init_err_msg, finl_err_msg
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH
      use geometry, only: imax2, jmax2, kmax2
      use geometry, only: imin1, jmin1, kmin1
      use geometry, only: imin2, kmin2

      implicit none

! variables related to global block structure
!---------------------------------------------------------------------//
      integer  :: dg_imin1, dg_imax1, &
                  dg_imin2, dg_imax2

      integer  :: dg_jmin1, dg_jmax1, &
                  dg_jmin2, dg_jmax2

      integer  :: dg_kmin1, dg_kmax1, &
                  dg_kmin2, dg_kmax2

! variables contain information of local indices of all processors
!---------------------------------------------------------------------//
      integer,dimension (:),allocatable ::      &
                  dg_istart1_all, dg_iend1_all, &
                  dg_istart2_all, dg_iend2_all, &
                  dg_isize_all

      integer,dimension (:),allocatable ::      &
                  dg_jstart1_all, dg_jend1_all, &
                  dg_jstart2_all, dg_jend2_all, &
                  dg_jsize_all

      integer,dimension (:),allocatable ::      &
                  dg_kstart1_all, dg_kend1_all, &
                  dg_kstart2_all, dg_kend2_all, &
                  dg_ksize_all

! variables relate to processor specific details
!---------------------------------------------------------------------//
      integer  :: dg_istart,  dg_iend,  &
                  dg_istart1, dg_iend1, &
                  dg_istart2, dg_iend2

      integer  :: dg_jstart,  dg_jend,  &
                  dg_jstart1, dg_jend1, &
                  dg_jstart2, dg_jend2

      integer  :: dg_kstart,  dg_kend,  &
                  dg_kstart1, dg_kend1, &
                  dg_kstart2, dg_kend2

      integer  :: dg_ijkstart2, dg_ijkend2, dg_ijksize2


! variables related to processor domain size
!---------------------------------------------------------------------//
      double precision :: dg_xstart, dg_xend, dg_dxinv
      double precision :: dg_ystart, dg_yend, dg_dyinv
      double precision :: dg_zstart, dg_zend, dg_dzinv

! Variables related to cyclic boundary  used in desgrid_functions
      integer,dimension(:,:),allocatable :: dg_cycoffset, icycoffset

      integer :: dg_pic_max_init = 25

! constants required for functions computing local and global ijkvalues
      integer dg_c1_gl, dg_c2_gl, dg_c3_gl  ! global
      integer dg_c1_lo, dg_c2_lo, dg_c3_lo  ! local

      integer, dimension(:), allocatable :: dg_c1_all
      integer, dimension(:), allocatable :: dg_c2_all
      integer, dimension(:), allocatable :: dg_c3_all

      contains

!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: PROCIJJK                                                  !
!                                                                      !
!  Purpose: Calculate the process that owns the cell with the given    !
!  I, J, K indicies.                                                   !
!......................................................................!
      integer function procijk(fi,fj,fk)
      implicit none
      integer fi,fj,fk
      procijk =fi+fj*nodesi+fk*nodesi*nodesj
      end function procijk


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: IofPROC                                                   !
!                                                                      !
!  Purpose: Return the I index of the current rank given an IJK value. !
!......................................................................!
      integer function iofproc(fijk)
      implicit none
      integer fijk
      iofproc = mod(mod(fijk,nodesi*nodesj),nodesi)
      end function iofproc


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: JofPROC                                                   !
!                                                                      !
!  Purpose: Return the J index of the current rank given an IJK value. !
!......................................................................!
      integer function jofproc(fijk)
      implicit none
      integer fijk
      jofproc = mod(fijk - iofproc(fijk),nodesi*nodesj)/nodesi
      end function jofproc


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: KofPROC                                                   !
!                                                                      !
!  Purpose: Return the K index of the current rank given an IJK value. !
!......................................................................!
      integer function kofproc(fijk)
      implicit none
      integer fijk
      kofproc = (fijk - iofproc(fijk) - jofproc(fijk)*nodesi) / &
         (nodesi*nodesj)
      end function kofproc

!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNIJK                                                 !
!                                                                      !
!  Purpose: Calculate the desgrid IJK given I, J, K components. The    !
!  result is local to the rank calculating the IJK.                    !
!......................................................................!
      integer function dg_funijk(fi,fj,fk)
      implicit none
      integer fi,fj,fk
      dg_funijk = fj+fi*dg_c1_lo+fk*dg_c2_lo+dg_c3_lo
      end function dg_funijk


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNIJK_GL                                              !
!                                                                      !
!  Purpose: Calculate the global desgrid IJK given I, J, K components. !
!  result should be consistent across all ranks.                       !
!......................................................................!
      integer function dg_funijk_gl(fi,fj,fk)
      implicit none
      integer fi,fj,fk
      dg_funijk_gl = fj+fi*dg_c1_gl+fk*dg_c2_gl+dg_c3_gl
      end function dg_funijk_gl


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNIJK_PROC                                            !
!                                                                      !
!  Purpose: Calculate the local desgrid IJK given I, J, K components   !
!  for the specified processor.                                        !
!......................................................................!
      integer function dg_funijk_proc(fi,fj,fk,fproc)
      implicit none
      integer fi,fj,fk,fproc
      dg_funijk_proc = fj + fi*dg_c1_all(fproc) + &
         fk*dg_c2_all(fproc) + dg_c3_all(fproc)
      end function dg_funijk_proc


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNIM                                                  !
!                                                                      !
!  Purpose: Return the local IJK value for the cell to the west.       !
!  This is equivalent to calculating: dg_funijk(I-1,J,K)               !
!......................................................................!
      integer function dg_funim(fijk)
      implicit none
      integer fijk
      dg_funim = fijk - dg_c1_lo
      end function dg_funim


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNIP                                                  !
!                                                                      !
!  Purpose: Return the local IJK value for the cell to the east.       !
!  This is equivalent to calculating: dg_funijk(I+1,J,K)               !
!......................................................................!
      integer function dg_funip(fijk)
      implicit none
      integer fijk
      dg_funip = fijk + dg_c1_lo
      end function dg_funip


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNJM                                                  !
!                                                                      !
!  Purpose: Return the local IJK value for the cell to the south.      !
!  This is equivalent to calculating: dg_funijk(I,J-1,K)               !
!......................................................................!
      integer function dg_funjm(fijk)
      implicit none
      integer fijk
      dg_funjm = fijk - 1
      end function dg_funjm


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNJP                                                  !
!                                                                      !
!  Purpose: Return the local IJK value for the cell to the north.      !
!  This is equivalent to calculating: dg_funijk(I,J+1,K)               !
!......................................................................!
      integer function dg_funjp(fijk)
      implicit none
      integer fijk
      dg_funjp = fijk + 1
      end function dg_funjp


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNKM                                                  !
!                                                                      !
!  Purpose: Return the local IJK value for the cell to the bottom.     !
!  This is equivalent to calculating: dg_funijk(I,J,K-1)               !
!......................................................................!
      integer function dg_funkm(fijk)
      implicit none
      integer fijk
      dg_funkm = fijk - dg_c2_lo
      end function dg_funkm


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_FUNJM                                                  !
!                                                                      !
!  Purpose: Return the local IJK value for the cell to the top.        !
!  This is equivalent to calculating: dg_funijk(I,J,K+1)               !
!......................................................................!
      integer function dg_funkp(fijk)
      implicit none
      integer fijk
      dg_funkp = fijk + dg_c2_lo
      end function dg_funkp


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_JOF_GL                                                 !
!                                                                      !
!  Purpose: Calculate the J index from a global IJK value.             !
!......................................................................!
      integer function dg_jof_gl(fijk)
      implicit none
      integer fijk
      dg_jof_gl = mod(mod(fijk-1,dg_c2_gl),dg_c1_gl)+dg_jmin2
      end function dg_jof_gl


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_IOF_GL                                                 !
!                                                                      !
!  Purpose: Calculate the I index from a global IJK value.             !
!......................................................................!
      integer function dg_iof_gl(fijk)
      implicit none
      integer fijk
      dg_iof_gl = (mod(fijk-dg_jof_gl(fijk)+dg_jmin2-1,dg_c2_gl)) /    &
         dg_c1_gl + dg_imin2
      end function dg_iof_gl


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_KOF_GL                                                 !
!                                                                      !
!  Purpose: Calculate the K index from a global IJK value.             !
!......................................................................!
      integer function dg_kof_gl(fijk)
      implicit none
      integer fijk
      dg_kof_gl = (fijk-dg_c3_gl-dg_iof_gl(fijk)*                      &
         dg_c1_gl-dg_jof_gl(fijk))/dg_c2_gl
      end function dg_kof_gl

!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_JOF_LO                                                 !
!                                                                      !
!  Purpose: Calculate the J index from a local IJK value.              !
!......................................................................!
      integer function dg_jof_lo(fijk)
      implicit none
      integer fijk
      dg_jof_lo = mod(mod(fijk-1,dg_c2_lo),dg_c1_lo)+dg_jstart2
      end function dg_jof_lo


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_IOF_LO                                                 !
!                                                                      !
!  Purpose: Calculate the I index from a local IJK value.              !
!......................................................................!
      integer function dg_iof_lo(fijk)
      implicit none
      integer fijk
      dg_iof_lo = (mod(fijk-dg_jof_lo(fijk)+dg_jstart2-1,dg_c2_lo)) /  &
         dg_c1_lo + dg_istart2
      end function dg_iof_lo


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_KOF_LO                                                 !
!                                                                      !
!  Purpose: Calculate the K index from a local IJK value.              !
!......................................................................!
      integer function dg_kof_lo(fijk)
      implicit none
      integer fijk
      dg_kof_lo = (fijk-dg_c3_lo-dg_iof_lo(fijk)*                      &
         dg_c1_lo-dg_jof_lo(fijk))/dg_c2_lo
      end function dg_kof_lo


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: DG_IJKCONV                                                !
!                                                                      !
!  Purpose: Convert the IJK from once process to another.              !
!                                                                      !
!......................................................................!
      integer function dg_ijkconv(fijk,fface,fto_proc)
      implicit none
      integer fijk,fto_proc,fface
      dg_ijkconv = dg_funijk_proc(dg_iof_lo(fijk)+                    &
        dg_cycoffset(fface,1), dg_jof_lo(fijk)+dg_cycoffset(fface,2), &
        dg_kof_lo(fijk)+dg_cycoffset(fface,3), fto_proc)
      end function dg_ijkconv


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: IofPOS                                                    !
!                                                                      !
!  Purpose: Calculate the cell I index containing the given position.  !
!......................................................................!
      integer function iofpos(fpos)
      implicit none
      double precision fpos
      iofpos = floor((fpos-dg_xstart)*dg_dxinv) + dg_istart1
      end function iofpos


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: JofPOS                                                    !
!                                                                      !
!  Purpose: Calculate the cell J index containing the given position.  !
!......................................................................!
      integer function jofpos(fpos)
      implicit none
      double precision fpos
      jofpos = floor((fpos-dg_ystart)*dg_dyinv) + dg_jstart1
      end function jofpos


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: KofPOS                                                    !
!                                                                      !
!  Purpose: Calculate the cell K index containing the given position.  !
!......................................................................!
      integer function kofpos(fpos)
      implicit none
      double precision fpos
      kofpos = floor((fpos-dg_zstart)*dg_dzinv) + dg_kstart1
      end function kofpos


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
!  Function: dg_is_ON_myPE_plus1layers                                 !
!                                                                      !
!  Purpose: Determine if the current rank contains the current cell    !
!  as its own or in its single layer of desgrid ghost cells.           !
!......................................................................!
      logical function dg_is_ON_myPE_plus1layers(lI, lJ, lK)
      implicit none
      integer, intent(in) :: lI, lJ, lK

      dg_is_ON_myPE_plus1layers = (&
         (dg_istart2 <= lI) .AND. (lI <= dg_iend2) .AND. &
         (dg_jstart2 <= lJ) .AND. (lJ <= dg_jend2) .AND. &
         (dg_kstart2 <= lK) .AND. (lK <= dg_kend2))

      end function


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DESGRID_INIT                                            !
!  Author: Pradeep G                                                   !
!                                                                      !
!  Purpose: Sets indices for the DESGRID and defines constants         !
!  required for the DESGRID functions. This is needed for              !
!  communication between the DESGRID for ghost cell exchanges.         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      subroutine desgrid_init()

      implicit none
!-----------------------------------------------
! Local varibles
!-----------------------------------------------
      double precision :: ltempdx,ltempdy,ltempdz
      integer :: lijkproc,liproc,ljproc,lkproc
      integer :: lijk
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("DESGRID_INIT")

! set indices for all processors
     allocate (dg_istart1_all(0:numpes-1), dg_iend1_all(0:numpes-1))
     allocate (dg_istart2_all(0:numpes-1), dg_iend2_all(0:numpes-1))
     allocate (dg_isize_all(0:nodesi-1))

     allocate (dg_jstart1_all(0:numpes-1), dg_jend1_all(0:numpes-1))
     allocate (dg_jstart2_all(0:numpes-1), dg_jend2_all(0:numpes-1))
     allocate (dg_jsize_all(0:nodesj-1))

     allocate (dg_kstart1_all(0:numpes-1), dg_kend1_all(0:numpes-1))
     allocate (dg_kstart2_all(0:numpes-1), dg_kend2_all(0:numpes-1))
     allocate (dg_ksize_all(0:nodesk-1))


      dg_istart1_all=0; dg_iend1_all=0
      dg_istart2_all=0; dg_iend2_all=0
      dg_isize_all=0

      dg_jstart1_all=0; dg_jend1_all=0
      dg_jstart2_all=0; dg_jend2_all=0
      dg_jsize_all=0

      dg_kstart1_all=0; dg_kend1_all=0
      dg_kstart2_all=0; dg_kend2_all=0
      dg_ksize_all=0


! set grid size based on user input desgridsearch_<ijk>max
      ltempdx = xlength/desgridsearch_imax
      ltempdy = ylength/desgridsearch_jmax
      ltempdz = zlength/desgridsearch_kmax

      dg_ksize_all(:) = 1

      lijkproc = 0
      do lkproc=0,nodesk-1
         do ljproc=0,nodesj-1
            do liproc=0,nodesi-1
               dg_isize_all(liproc) = NINT((xe(iend1_all(lijkproc))-xe(istart1_all(lijkproc)-1))/ltempdx)
               dg_jsize_all(ljproc) = NINT((yn(jend1_all(lijkproc))-yn(jstart1_all(lijkproc)-1))/ltempdy)
               dg_ksize_all(lkproc) = NINT((zt(kend1_all(lijkproc))-zt(kstart1_all(lijkproc)-1))/ltempdz)
               dg_istart1_all(lijkproc) = sum(dg_isize_all(0:liproc-1)) + 2
               dg_jstart1_all(lijkproc) = sum(dg_jsize_all(0:ljproc-1)) + 2
               dg_kstart1_all(lijkproc) = sum(dg_ksize_all(0:lkproc-1)) + 2
               dg_iend1_all(lijkproc) = dg_isize_all(liproc)+dg_istart1_all(lijkproc)-1
               dg_jend1_all(lijkproc) = dg_jsize_all(ljproc)+dg_jstart1_all(lijkproc)-1
               dg_kend1_all(lijkproc) = dg_ksize_all(lkproc)+dg_kstart1_all(lijkproc)-1
               dg_istart2_all(lijkproc) = dg_istart1_all(lijkproc)-1
               dg_jstart2_all(lijkproc) = dg_jstart1_all(lijkproc)-1
               dg_kstart2_all(lijkproc) = dg_kstart1_all(lijkproc)-1
               dg_iend2_all(lijkproc) =  dg_iend1_all(lijkproc)+1
               dg_jend2_all(lijkproc) = dg_jend1_all(lijkproc)+1
               dg_kend2_all(lijkproc) = dg_kend1_all(lijkproc)+1
               lijkproc = lijkproc + 1
            end do
         end do
      end do

! set indices for global block
      dg_imin2 = 1
      dg_imin1 = 2
      dg_imax1 = dg_imin1+sum(dg_isize_all(0:nodesi-1))-1
      dg_imax2 = dg_imax1+1
      dg_jmin2 = 1
      dg_jmin1 = 2
      dg_jmax1 = dg_jmin1+sum(dg_jsize_all(0:nodesj-1))-1
      dg_jmax2 = dg_jmax1+1
      dg_kmin2 = 1
      dg_kmin1 = 2
      dg_kmax1 = dg_kmin1+sum(dg_ksize_all(0:nodesk-1))-1
      dg_kmax2 = dg_kmax1+1

! set offset indices for periodic boundaries
      lijkproc = mype
      liproc = iofproc(lijkproc)
      ljproc = jofproc(lijkproc)
      lkproc = kofproc(lijkproc)
      allocate(dg_cycoffset(2*dimn,3),icycoffset(2*dimn,3))
      dg_cycoffset(:,:) = 0; icycoffset(:,:) = 0
      if (des_periodic_walls_x) then
         if(liproc.eq.0) then
            dg_cycoffset(2,1)= (dg_imax2-dg_imin1)
            icycoffset(2,1)= (imax2-imin1)
         end if
         if(liproc.eq.nodesi-1) then
            dg_cycoffset(1,1)=-(dg_imax2-dg_imin1)
            icycoffset(1,1)= -(imax2-imin1)
         end if
      end if
      if (des_periodic_walls_y) then
         if(ljproc.eq.0) then
            dg_cycoffset(4,2)= (dg_jmax2-dg_jmin1)
            icycoffset(4,2)= (jmax2-jmin1)
         end if
         if(ljproc.eq.nodesj-1) then
            dg_cycoffset(3,2)=-(dg_jmax2-dg_jmin1)
            icycoffset(3,2)= -(jmax2-jmin1)
         end if
      end if
      if (des_periodic_walls_z) then
         if(lkproc.eq.0) then
            dg_cycoffset(6,3)=(dg_kmax2-dg_kmin1)
            icycoffset(6,3)= (kmax2-kmin1)
         end if
         if(lkproc.eq.nodesk-1) then
            dg_cycoffset(5,3)=-(dg_kmax2-dg_kmin1)
            icycoffset(5,3)= -(kmax2-kmin1)
         end if
      end if


! set the indices and variables for the current processor
      dg_istart2 = dg_istart2_all(mype)
      dg_istart1 = dg_istart1_all(mype)
      dg_iend1 = dg_iend1_all(mype)
      dg_iend2 = dg_iend2_all(mype)
      dg_jstart2 = dg_jstart2_all(mype)
      dg_jstart1 = dg_jstart1_all(mype)
      dg_jend1 = dg_jend1_all(mype)
      dg_jend2 = dg_jend2_all(mype)
      dg_kstart2 = dg_kstart2_all(mype)
      dg_kstart1 = dg_kstart1_all(mype)
      dg_kend1 = dg_kend1_all(mype)
      dg_kend2 = dg_kend2_all(mype)

!       Defining new set of varaibles to define upper and lower bound of the
! indices to include actual physical boundaries of the problem.
      dg_istart = dg_istart1;   dg_iend = dg_iend1
      dg_jstart = dg_jstart1;   dg_jend = dg_jend1
      dg_kstart = dg_kstart1;   dg_kend = dg_kend1

      if(dg_istart .eq. dg_imin1) dg_istart = dg_imin2
      if(dg_iend   .eq. dg_imax1) dg_iend   = dg_imax2

      if(dg_jstart .eq. dg_jmin1) dg_jstart = dg_jmin2
      if(dg_jend   .eq. dg_jmax1) dg_jend   = dg_jmax2

      if(dg_kstart .eq. dg_kmin1) dg_kstart = dg_kmin2
      if(dg_kend   .eq. dg_kmax1) dg_kend   = dg_kmax2

! set constants required for dg_funijk dg_funijk_gl for all processors
      allocate(dg_c1_all(0:numpes-1),dg_c2_all(0:numpes-1),dg_c3_all(0:numpes-1))

      dg_c1_all=0;dg_c2_all=0;dg_c3_all=0
      do lijkproc = 0,numpes-1
         dg_c1_all(lijkproc) = (dg_jend2_all(lijkproc)-dg_jstart2_all(lijkproc)+1)
         dg_c2_all(lijkproc) = dg_c1_all(lijkproc)*(dg_iend2_all(lijkproc)-dg_istart2_all(lijkproc)+1)
         dg_c3_all(lijkproc) = -dg_c1_all(lijkproc)*dg_istart2_all(lijkproc) &
                                  -dg_c2_all(lijkproc)*dg_kstart2_all(lijkproc) &
                                  -dg_jstart2_all(lijkproc)+1
      end do
! set global constants
      dg_c1_gl = (dg_jmax2-dg_jmin2+1)
      dg_c2_gl = dg_c1_gl*(dg_imax2-dg_imin2+1)
      dg_c3_gl = -dg_c1_gl*imin2-dg_c2_gl*kmin2-dg_jmin2+1

! set local constants
      dg_c1_lo = dg_c1_all(mype)
      dg_c2_lo = dg_c2_all(mype)
      dg_c3_lo = dg_c3_all(mype)

      dg_ijksize2 = (dg_iend2-dg_istart2+1)* &
                    (dg_jend2-dg_jstart2+1)* &
                    (dg_kend2-dg_kstart2+1)
      dg_ijkstart2 = dg_funijk(dg_istart2,dg_jstart2,dg_kstart2)
      dg_ijkend2 = dg_funijk(dg_iend2,dg_jend2,dg_kend2)

! Confirmation checks
      IF (DG_IJKSTART2.NE.1) THEN
         WRITE(ERR_MSG,1100)'DG_IJKStart2'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(DG_IJKEND2 /= DG_IJKSIZE2) THEN
         WRITE(ERR_MSG,1100)'DG_IJKEnd2'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Invalid DG_IJKStart2. FATAL')

      IF(DG_IMIN1 > DG_IMAX1) THEN
         WRITE(ERR_MSG,1105) 'DG_IMIN1 > DG_IMAX1'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(DG_JMIN1 > DG_JMAX1) THEN
         WRITE(ERR_MSG,1105) 'DG_JMIN1 > DG_JMAX1'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(DG_KMIN1 > DG_KMAX1) THEN
         WRITE(ERR_MSG,1105) 'DG_KMIN1 > DG_KMAX1'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1105 FORMAT('Error 1105: Invalid DES grid indices: ',A,/'This is ',   &
         'likely the result of automated grid calculations based',/    &
         'on the maximum particle size. Specify the number of ',       &
         'partitions',/'for the DES grid in the mfix.dat file ',       &
         '(e.g., DESGRIDSEARCH_IMAX)')


! set the domain length and dx,dy and dz values used in particle_in_cell
! to bin the particles
      dg_xstart = xe(istart2)
      dg_xend = xe(iend1)
      dg_ystart = yn(jstart2)
      dg_yend = yn(jend1)
      dg_dxinv = (dg_iend1-dg_istart1+1)/(dg_xend-dg_xstart)
      dg_dyinv = (dg_jend1-dg_jstart1+1)/(dg_yend-dg_ystart)
      dg_zstart = zt(kstart2)
      dg_zend = zt(kend1)
      dg_dzinv = (dg_kend1-dg_kstart1+1)/(dg_zend-dg_zstart)

! allocate the desgridsearch related variables
      allocate(dg_pic(dg_ijksize2))
      dg_pic(:)%isize = 0
      do lijk = 1,dg_ijksize2
         allocate(dg_pic(lijk)%p(dg_pic_max_init))
      end do

!      call des_dbggrid
      CALL FINL_ERR_MSG
      end subroutine desgrid_init

!------------------------------------------------------------------------
! Subroutine       : desgrid_pic
! Purpose          : it updates the particle in cell information
!                    and also located the particle
!                    moving across processor boundary
! Parameters       : plocate - Locate the particles
! Comments         : plocate should be set to true only when particle has to
!                    be located; if it is false only PIC information will be
!                    updated by this routine
!------------------------------------------------------------------------
      subroutine desgrid_pic(plocate, dg_pijkprv, dg_pijk, particle_state, des_pos_new)

      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      logical, INTENT(IN) :: plocate

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_pos_new
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_state
      INTEGER, DIMENSION(:), INTENT(INOUT) :: dg_pijk
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijkprv

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer, dimension(dg_ijksize2) :: lpic,lindx
      integer li,lj,lk,lijk,lijk_count,lcurpar,lparcount,lcurpic
      logical, save :: first_pass = .true.
!-----------------------------------------------

! locate the particles including ghost cells
      max_isize = 0
      lparcount = 1
      lpic(:) = 0
      if (plocate) then
         dg_pijkprv(:)= dg_pijk(:)
         do lcurpar = 1,max_pip
            if(lparcount.gt.pip) exit
            if(nonexistent==particle_state(lcurpar)) cycle

            lparcount = lparcount + 1
            li = min(dg_iend2,max(dg_istart2,iofpos(des_pos_new(lcurpar,1))))
            lj = min(dg_jend2,max(dg_jstart2,jofpos(des_pos_new(lcurpar,2))))
            lk = min(dg_kend2,max(dg_kstart2,kofpos(des_pos_new(lcurpar,3))))
            dg_pijk(lcurpar) = dg_funijk(li,lj,lk)
            lijk = dg_pijk(lcurpar)
            lpic(lijk) = lpic(lijk) + 1
         end do
      else
         do lcurpar = 1,max_pip
            if(lparcount.gt.pip) exit
            if(nonexistent==particle_state(lcurpar)) cycle
            lparcount = lparcount + 1
            lijk = dg_pijk(lcurpar)
            lpic(lijk) = lpic(lijk) + 1
         end do
      end if

      if (first_pass) then
         dg_pijkprv(:) = dg_pijk(:)
         first_pass = .false.
      end if

! redefine the array of dg_pic
      do lijk = dg_ijkstart2,dg_ijkend2
         lcurpic = lpic(lijk)
         if(lcurpic > size(dg_pic(lijk)%p)) then
            deallocate(dg_pic(lijk)%p)
            allocate(dg_pic(lijk)%p(2*lcurpic))
         end if
         dg_pic(lijk)%isize = lcurpic
         max_isize = max(max_isize,dg_pic(lijk)%isize)
      end do

! assign the particle info in pic array
      lindx(:) = 1
      lparcount = 1

      do lcurpar = 1, max_pip
         if(lparcount.gt.pip) cycle
         if(nonexistent==particle_state(lcurpar)) cycle
         lijk = dg_pijk(lcurpar)

         lijk_count = lindx(lijk)
         lindx(lijk) = lindx(lijk) +  1

         dg_pic(lijk)%p(lijk_count) = lcurpar

         lparcount = lparcount + 1
      enddo


!      open (unit=100,file='desgrid.txt',status='unknown')
!      write(100,*)lpic
!      close(100)

    contains

      include 'functions.inc'

    end subroutine desgrid_pic

!------------------------------------------------------------------------
! subroutine       : desgrid_neigh_build ()
! Purpose          : This particles build the neigh list for the particles
!                    currently active in the system
!------------------------------------------------------------------------
      subroutine desgrid_neigh_build(des_pos_new, dg_pijk, particle_state, des_radius, neighbor_index)

      implicit none

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: des_radius
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_pos_new
      INTEGER, DIMENSION(:), INTENT(INOUT) :: NEIGHBOR_INDEX
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_state
      INTEGER, DIMENSION(:), INTENT(IN) :: dg_pijk

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lcurpar,lkoffset
      integer lijk,lic,ljc,lkc,li,lj,lk,ltotpic,lpicloc,lneigh,cc
      double precision lsearch_rad,ldistsquared
      double precision :: ldistvec(3)
      double precision :: lcurpar_pos(3)
      double precision :: lcur_off
      integer il_off,iu_off,jl_off,ju_off,kl_off,ku_off
      integer, dimension(:), allocatable :: tmp_neigh

!-----------------------------------------------

! loop through neighbours and build the contact particles list for particles
! present in the system
      lkoffset = dimn-2

      allocate(tmp_neigh(max_isize))

      do lcurpar =1,max_pip

!$  if (.false.) then
            if (NEIGHBOR_INDEX(lcurpar) .eq. 0) then
                if (lcurpar .eq. 1) then
                    NEIGHBOR_INDEX(lcurpar) = 1
                else
                    NEIGHBOR_INDEX(lcurpar) = NEIGHBOR_INDEX(lcurpar-1)
                endif
            endif
!$  endif

            if (nonexistent==particle_state(lcurpar) .or.entering_particle==particle_state(lcurpar) &
                 .or. entering_ghost==particle_state(lcurpar) &
                 .or. normal_ghost==particle_state(lcurpar) &
                 .or. exiting_ghost==particle_state(lcurpar)) cycle
         lijk = dg_pijk(lcurpar)
         lic = dg_iof_lo(lijk)
         ljc = dg_jof_lo(lijk)
         lkc = dg_kof_lo(lijk)

         il_off = 1
         iu_off = 1
         jl_off = 1
         ju_off = 1
         kl_off = 1
         ku_off = 1

         lcurpar_pos(:) = des_pos_new(lcurpar,:)
!   The desgrid size should not be less than 2*dia*rlm_factor
         lcur_off = (lcurpar_pos(1)-dg_xstart)*dg_dxinv - &
            floor((lcurpar_pos(1)-dg_xstart)*dg_dxinv)
         if(lcur_off .ge. 0.5) then
            il_off = 0
         else
            iu_off = 0
         endif

         lcur_off = (lcurpar_pos(2)-dg_ystart)*dg_dyinv - &
            floor((lcurpar_pos(2)-dg_ystart)*dg_dyinv)
         if(lcur_off .ge. 0.5) then
            jl_off = 0
         else
            ju_off = 0
         endif

         lcur_off = (lcurpar_pos(3)-dg_zstart)*dg_dzinv - &
            floor((lcurpar_pos(3)-dg_zstart)*dg_dzinv)
         if(lcur_off .ge. 0.5) then
            kl_off = 0
         else
            ku_off = 0
         endif

         do lk = lkc-kl_off,lkc+ku_off
         do lj = ljc-jl_off,ljc+ju_off
         do li = lic-il_off,lic+iu_off
            lijk = dg_funijk(li,lj,lk)
            ltotpic =dg_pic(lijk)%isize
            do lpicloc = 1,ltotpic
               lneigh = dg_pic(lijk)%p(lpicloc)

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

    contains

      include 'functions.inc'

    end subroutine desgrid_neigh_build

!------------------------------------------------------------------------
! subroutine       : des_dbggrid
! Purpose          : For printing the indices used for desgrid
!------------------------------------------------------------------------
      subroutine des_dbggrid()

      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lproc,liproc,ljproc,lkproc
      character (255) filename
!-----------------------------------------------

      write(filename,'("dbg_desgridn",I4.4,".dat")') mype
      open(44,file=filename)
      do lproc = 0,numpes-1
         write(44,*) "Information for Proc =", lproc
         liproc= iofproc(lproc)
         ljproc= jofproc(lproc)
         lkproc= kofproc(lproc)
         write(44,*) "   "
         write(44,*) "i,j,k location of proc",liproc,ljproc,lkproc
         write(44,*) "i,j,k size of proc", dg_isize_all(liproc),dg_jsize_all(ljproc),dg_ksize_all(lkproc)
         write(44,*) "-------------------------------------------------"
         write(44,*) "indices        start     end"
         write(44,*) "-------------------------------------------------"
         write(44,*) "for i1:      ",dg_istart1_all(lproc),dg_iend1_all(lproc)
         write(44,*) "for i2:      ",dg_istart2_all(lproc),dg_iend2_all(lproc)
         write(44,*) "for j1:      ",dg_jstart1_all(lproc),dg_jend1_all(lproc)
         write(44,*) "for j2:      ",dg_jstart2_all(lproc),dg_jend2_all(lproc)
         write(44,*) "for k1:      ",dg_kstart1_all(lproc),dg_kend1_all(lproc)
         write(44,*) "for k2:      ",dg_kstart2_all(lproc),dg_kend2_all(lproc)
      end do
      write(44,*) "   "
      write(44,*) "-------------------------------------------------"
      write(44,*) "Local Start and end"
      write(44,*) "-------------------------------------------------"
      write(44,*) "for i :      ",dg_istart, dg_iend
      write(44,*) "for i1:      ",dg_istart1,dg_iend1
      write(44,*) "for i2:      ",dg_istart2,dg_iend2
      write(44,*) "   "
      write(44,*) "for j :      ",dg_jstart, dg_jend
      write(44,*) "for j1:      ",dg_jstart1,dg_jend1
      write(44,*) "for j2:      ",dg_jstart2,dg_jend2
      write(44,*) "   "
      write(44,*) "for k :      ",dg_kstart, dg_kend
      write(44,*) "for k1:      ",dg_kstart1,dg_kend1
      write(44,*) "for k2:      ",dg_kstart2,dg_kend2
      write(44,*) "   "
      write(44,*) "-------------------------------------------------"
      write(44,*) "global Start and end"
      write(44,*) "-------------------------------------------------"
      write(44,*) "for i1:      ",dg_imin1,dg_imax1
      write(44,*) "for i2:      ",dg_imin2,dg_imax2
      write(44,*) "for j1:      ",dg_jmin1,dg_jmax1
      write(44,*) "for j2:      ",dg_jmin2,dg_jmax2
      write(44,*) "for k1:      ",dg_kmin1,dg_kmax1
      write(44,*) "for k2:      ",dg_kmin2,dg_kmax2
      write(44,*) "   "
      write(44,*) "-------------------------------------------------"
      write(44,*) "dg_xstart:   ",dg_xstart
      write(44,*) "dg_xend:     ",dg_xend
      write(44,*) "dg_dxinv:    ",dg_dxinv
      write(44,*) "   "
      write(44,*) "dg_ystart:   ",dg_ystart
      write(44,*) "dg_yend:     ",dg_yend
      write(44,*) "dg_dyinv:    ",dg_dyinv
      write(44,*) "   "
      write(44,*) "dg_zstart:   ",dg_zstart
      write(44,*) "dg_zend:     ",dg_zend
      write(44,*) "dg_dzinv:    ",dg_dzinv
      write(44,*) "   "
      write(44,*) "dg_c1_lo/gl: ",dg_c1_lo, dg_c1_gl
      write(44,*) "dg_c2_lo/gl: ",dg_c2_lo, dg_c2_gl
      write(44,*) "dg_c3_lo/gl: ",dg_c3_lo, dg_c3_gl
      write(44,*) "   "
      write(44,*) "-------------------------------------------------"

      close(44)
      end subroutine des_dbggrid

      end module
