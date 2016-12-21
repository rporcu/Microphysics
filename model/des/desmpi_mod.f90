!----------------------------------------------------------------------!
!  Module: DESMPI                                                      !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for packing real and ghost particles     !
!     into the MPI send buffers.                                       !
!----------------------------------------------------------------------!
      module desmpi

! Ghost particle packet size.
      INTEGER :: iGhostPacketSize
      INTEGER :: iParticlePacketSize
      INTEGER :: iPairPacketSize

! Flags and constants for interfaces
      integer, dimension(:), allocatable :: ineighproc
      logical, dimension(:), allocatable :: iexchflag

! offset for periodic boundaries
      double precision, dimension(:,:), allocatable :: dcycl_offset

      type array
         double precision, dimension(:), allocatable :: facebuf
      end type array

! following variables used for sendrecv ghost particles and particle exchange
      type(array), dimension(:), allocatable :: dsendbuf
      type(array), dimension(:), allocatable :: drecvbuf

      integer,dimension(:),allocatable:: isendcnt
      integer,dimension(:),allocatable:: isendreq
      integer,dimension(:),allocatable:: irecvreq

      integer,parameter :: ibufoffset = 2

! The maximum size of the receive buffer.
      integer :: imaxbuf
      integer :: ispot

! following variables are used for gather and scatter
      double precision, dimension(:), allocatable :: drootbuf
      double precision, dimension(:), allocatable :: dprocbuf
      integer, dimension(:), allocatable :: irootbuf
      integer, dimension(:), allocatable :: iprocbuf

      integer,dimension(:), allocatable:: idispls
      integer,dimension(:), allocatable:: iscattercnts
      integer,dimension(:), allocatable:: igathercnts

      integer :: iscr_recvcnt
      integer :: igath_sendcnt

! following variables are used to identify the cell number for ghost cells
      integer,dimension(:,:),allocatable :: isendindices
      integer,dimension(:,:),allocatable :: irecvindices

! variables used to read initial particle properties
      double precision, dimension(:,:), allocatable:: dpar_pos
      double precision, dimension(:,:), allocatable:: dpar_vel
      double precision, dimension(:), allocatable:: dpar_den
      double precision, dimension(:), allocatable:: dpar_rad

      contains

!------------------------------------------------------------------------
! subroutine       : des_dbgmpi
! Purpose          : For printing the flags and values set for interface
!                    communication
! Parameters       : ptype - based on this following info is printed to
!                    the file
!                    1 - interface flags
!                    2 - send buffer for ghost particles
!                    3 - recv buffer for ghost particles
!                    4 - particle information
!                    5 - send buffer for particles exchanging processor
!                    6 - particles info
!                    7 - neighinfo
!------------------------------------------------------------------------
      subroutine des_dbgmpi(ptype)
      implicit none
      integer ptype

      end subroutine des_dbgmpi


      end module
