!----------------------------------------------------------------------!
!  Module: MPI_INIT_DES                                                !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for setting up DES MPI communications.   !
!                                                                      !
!----------------------------------------------------------------------!
      module mpi_init_des

      use compar, only: iend1_all, jend1_all, kend1_all
      use compar, only: istart1_all, jstart1_all, kstart1_all
      use compar, only: mype, pe_io, numpes
      use compar, only: nodesi, nodesj, nodesk
      use desgrid, only: desgrid_pic, iofproc, jofproc, kofproc, procijk, dg_funijk
      use desgrid, only: dg_iend, dg_jend, dg_kend
      use desgrid, only: dg_iend1, dg_jend1, dg_kend1
      use desgrid, only: dg_iend2, dg_jend2, dg_kend2
      use desgrid, only: dg_imax1, dg_jmax1, dg_kmax1
      use desgrid, only: dg_imin1, dg_jmin1, dg_kmin1
      use desgrid, only: dg_istart, dg_jstart, dg_kstart
      use desgrid, only: dg_istart1, dg_jstart1, dg_kstart1
      use desgrid, only: dg_istart2, dg_jstart2, dg_kstart2
      use desmpi, only: dcycl_offset, ispot, imaxbuf, dpar_rad, dpar_vel, dpar_den
      use desmpi, only: dprocbuf, iprocbuf, drootbuf, irootbuf, dpar_pos, idispls, iscr_recvcnt, iscattercnts, ineighproc
      use desmpi, only: drecvbuf, dsendbuf, isendcnt, isendindices, irecvindices, isendreq, irecvreq, iexchflag, igathercnts
      use discretelement, only: des_periodic_walls_x, des_periodic_walls_y, des_periodic_walls_z
      use discretelement, only: normal_particle, do_nsearch, do_old, des_explicitly_coupled
      use discretelement, only: particles, dimn, ighost_updated, pip, max_pip, xe, yn, zt
      use error_manager, only: err_msg, init_err_msg, flush_err_msg, finl_err_msg
      use param, only: DIMENSION_N_s

      contains

!----------------------------------------------------------------------!
!  Module: DESMPI_INIT                                                 !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Allocates and initializes variables used by MPI send/recv  !
!  calls. Sets flags related to periodic boundaries and processor      !
!  interfaces.                                                         !
!                                                                      !
!----------------------------------------------------------------------!
      subroutine desmpi_init

      use desmpi, only: iGhostPacketSize
      use desmpi, only: iParticlePacketSize
      use desmpi, only: iPairPacketSize
      use discretelement, only: DES_USR_VAR_SIZE

      implicit none

!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer :: lfaces
      integer :: lmaxlen1,lmaxlen2,lmaxarea,lmaxghostpar,ii

      DOUBLE PRECISION, PARAMETER :: ONEMBo8 = 131072.0

!-----------------------------------------------

! Calculate the size of ghost particle packets:
      iGhostPacketSize = 15 + DES_USR_VAR_SIZE

! Calculate the size of particle packets.
      iParticlePacketSize = 30 + DES_USR_VAR_SIZE
      IF(DES_EXPLICITLY_COUPLED) &
         iParticlePacketSize = iParticlePacketSize + 3
      IF(DO_OLD) &
         iParticlePacketSize = iParticlePacketSize + 15

! Calculate the size of neighbor data
      iPairPacketSize = 11

! Calculate the initial size of send and recv buffer based on max_pip,
! total cells max number of boundary cells, and ghost par packet size.
      lfaces = dimn*2
      lmaxlen1 = dg_iend2-dg_istart2+1
      lmaxlen2 = dg_jend2-dg_jstart2+1
      lmaxlen1 = max(lmaxlen1,dg_kend2 -dg_kstart2+1)
      lmaxlen2 = max(lmaxlen2,dg_kend2 -dg_kstart2+1)


! Note: 10 is added for buffer and is required for send and recv indices
      lmaxarea = lmaxlen1*lmaxlen2 + 10

! Random value. This gets resized by DESMPI_CHECK_SENDRECVBUF
      lmaxghostpar = 100
      imaxbuf = lmaxghostpar*lmaxarea*iGhostPacketSize

      WRITE(ERR_MSG, 1000) iMAXBUF/ONEMBo8, ONEMBo8/iGhostPacketSize,  &
         ONEMBo8/iParticlePacketSize, ONEMBo8/iPairPacketSize
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT(/'DES MPI send/recv buffer: ',F7.1,' MB',/' o  ',F6.0,1X, &
         'Ghost Particles/MB',/' o  ',F6.0,1X,'Particles/MB',/' o  ',  &
         F6.0,1X,'Neighbor Pairs/MB')


      allocate (dsendbuf(2));
      allocate (drecvbuf(2));
      do ii=1, size(dsendbuf)
         allocate (dsendbuf(ii)%facebuf(imaxbuf));
         allocate (drecvbuf(ii)%facebuf(imaxbuf));
      end do

      allocate (isendindices(lmaxarea,lfaces)); isendindices=0
      allocate (irecvindices(lmaxarea,lfaces)); irecvindices=0

      allocate (isendreq(lfaces)); isendreq=0
      allocate (irecvreq(lfaces)); irecvreq=0
      allocate (isendcnt(lfaces)); isendcnt=0

      allocate (dcycl_offset(lfaces,dimn)); dcycl_offset=0.0
      allocate (ineighproc(lfaces)); ineighproc=0
      allocate (iexchflag(lfaces)); iexchflag=.FALSE.

! allocate variables related to scattergather
      allocate(iscattercnts(0:numpes-1)); iscattercnts=0
      allocate(igathercnts(0:numpes-1));  igathercnts=0
      allocate(idispls(0:numpes-1)); idispls=0

      end subroutine desmpi_init


      end module mpi_init_des
