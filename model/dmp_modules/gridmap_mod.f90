!--------------------------------------------------------------------
! Purpose:
! Contains following subroutines:
!    partition, gridmap_init
!--------------------------------------------------------------------
       module gridmap

!-----------------------------------------------
! Modules
!-----------------------------------------------
          use geometry, only: imax1, imin1, jmax1, jmin1, kmax1, kmin1
          use compar, only: nodesi, nodesj, nodesk
          use compar, only: nlayers_bicgs
          use compar, only: istart1_all, jstart1_all, kstart1_all, iend1_all, jend1_all, kend1_all
          use error_manager, only: err_msg, init_err_msg, flush_err_msg, finl_err_msg

        implicit none

        contains

!----------------------------------------------------------------------!
! Purpose: Routine to partition the grid. It works for 1-d, 2-d        !
! decomposition in the current implementation                          !
!----------------------------------------------------------------------!
      SUBROUTINE PARTITION(CYC_XLL, CYC_YLL, CYC_ZLL)

      implicit none

! DUMMY Arguments
!---------------------------------------------------------------------//
! Local flags for cyclic boundarys
       LOGICAL :: CYC_XLL, CYC_YLL, CYC_ZLL

! Local variables
!---------------------------------------------------------------------//
      INTEGER, dimension(0:nodesi-1) :: isize1_all
      INTEGER, dimension(0:nodesj-1) :: jsize1_all
      INTEGER, dimension(0:nodesk-1) :: ksize1_all

      INTEGER :: ip, iproc, isize, iremain
      INTEGER :: jp, jproc, jsize, jremain
      INTEGER :: kp, kproc, ksize, kremain

      INTEGER :: ijkproc

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("PARTITION")

! Set the number of layers for BICGSTAB
      IF(NODESI .NE. 1 .AND. CYC_XLL) nlayers_bicgs = 2
      IF(NODESJ .NE. 1 .AND. CYC_YLL) nlayers_bicgs = 2
      IF(NODESK .NE. 1 .AND. CYC_ZLL) nlayers_bicgs = 2

! Flag that the current setup may not be efficient.
      IF(NODESJ .NE. 1) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG
      ENDIF

! Determine the size in i direction and add the remainder sequentially
      isize = (imax1-imin1+1)/nodesi
      isize1_all(0:nodesi-1) = isize
      iremain = (imax1-imin1+1) - nodesi*isize
      IF (iremain.ge.1) isize1_all( 0:(iremain-1) ) = isize + 1

! Determine the size in j direction and add the remainder sequentially
      jsize = (jmax1-jmin1+1)/nodesj
      jsize1_all(0:nodesj-1) = jsize
      jremain = (jmax1-jmin1+1) - nodesj*jsize
      IF (jremain.ge.1) jsize1_all( 0:(jremain-1) ) = jsize + 1

! Determine the size in k direction and add the remainder sequentially
      ksize = (kmax1-kmin1+1)/nodesk
      ksize1_all(0:nodesk-1) = ksize
      kremain = (kmax1-kmin1+1) - nodesk*ksize
      IF (kremain.ge.1) ksize1_all( 0:(kremain-1) ) = ksize + 1


! The following is general for 1-d or 2-d or 3-d decompostion
! Determining  istart, jstart and kstart for all the processors
      ijkproc = 0
      kp = kmin1
      do kproc=0,nodesk-1
         jp = jmin1
         do jproc=0,nodesj-1
            ip = imin1
            do iproc=0,nodesi-1

               istart1_all(ijkproc) = ip + sum(isize1_all(0:iproc-1))
               iend1_all(ijkproc) = istart1_all(ijkproc) + isize1_all(iproc)-1

               jstart1_all(ijkproc) = jp + sum(jsize1_all(0:jproc-1))
               jend1_all(ijkproc) = jstart1_all(ijkproc) + jsize1_all(jproc)-1

               kstart1_all(ijkproc) = kp + sum(ksize1_all(0:kproc-1))
               kend1_all(ijkproc) = kstart1_all(ijkproc) + ksize1_all(kproc)-1

               ijkproc = ijkproc+1

            ENDDO
         ENDDO
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('WARNING 1000: The preconditioner for the linear solver',/&
         'MIGHT NOT be very efficient with DMP partitions in the y-',  &
         'axis.')

      END SUBROUTINE PARTITION


!----------------------------------------------------------------------!
! Purpose: Initializing all the variables from the information         !
! obtained in the above routine.                                       !
!----------------------------------------------------------------------!
        SUBROUTINE GRIDMAP_INIT

        use functions, only: funijk
        use compar, only: ijksize3_all, ijkstart3_all, ijkend3_all
        use compar, only: istart_all, iend_all
        use compar, only: jstart_all, jend_all
        use compar, only: kstart_all, kend_all
        use compar, only: istart2_all, iend2_all
        use compar, only: jstart2_all, jend2_all
        use compar, only: kstart2_all, kend2_all
        use compar, only: istart3_all, iend3_all
        use compar, only: jstart3_all, jend3_all
        use compar, only: kstart3_all, kend3_all
        use compar, only: imap, jmap, kmap
        use compar, only: imap_c, jmap_c, kmap_c
        use compar, only: imap_c, jmap_c, kmap_c
        use compar, only: c0, c1, c2
        use compar
        use geometry, only: imin2, jmin2, kmin2
        use geometry, only: imin3, jmin3, kmin3
        use geometry, only: imin4, jmin4, kmin4
        use geometry, only: imax2, jmax2, kmax2
        use geometry, only: imax3, jmax3, kmax3
        use geometry, only: imax4, jmax4, kmax4
        use geometry, only: cyclic_x, cyclic_y, cyclic_z
        use geometry, only: cyclic_x_pd, cyclic_y_pd, cyclic_z_pd

        implicit none

! Local variables
!---------------------------------------------------------------------//
! Loop indicies
      integer :: iproc, ii, jj, kk
! Local flags for cyclic boundarys
      LOGICAL :: CYC_XL, CYC_YL, CYC_ZL
! Amount of load imbalance
      INTEGER :: IMBALANCE
! Theoritical speedup (based on load imbalance)
      CHARACTER(len=32) :: AMDAHL_SPEEDUP

!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("GRIDMAP_INIT")

! Set local flags for cyclic boundaries.
      CYC_XL = (CYCLIC_X .OR. CYCLIC_X_PD)
      CYC_YL = (CYCLIC_Y .OR. CYCLIC_Y_PD)
      CYC_ZL = (CYCLIC_Z .OR. CYCLIC_Z_PD)

      IF(.NOT.ALLOCATED(ijksize3_all))   allocate( ijksize3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(ijkstart3_all))  allocate( ijkstart3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(ijkend3_all))    allocate( ijkend3_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart_all))     allocate( istart_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart_all))     allocate( jstart_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart_all))     allocate( kstart_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart1_all))    allocate( istart1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart1_all))    allocate( jstart1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart1_all))    allocate( kstart1_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart2_all))    allocate( istart2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart2_all))    allocate( jstart2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart2_all))    allocate( kstart2_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart3_all))    allocate( istart3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart3_all))    allocate( jstart3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart3_all))    allocate( kstart3_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend_all))       allocate( iend_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend_all))       allocate( jend_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend_all))       allocate( kend_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend1_all))      allocate( iend1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend1_all))      allocate( jend1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend1_all))      allocate( kend1_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend2_all))      allocate( iend2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend2_all))      allocate( jend2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend2_all))      allocate( kend2_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend3_all))      allocate( iend3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend3_all))      allocate( jend3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend3_all))      allocate( kend3_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(displs))         allocate( displs(0:numPEs-1) )


      CALL PARTITION(CYC_XL, CYC_YL, CYC_ZL)

! The upper and lower bounds are prescribed such that two ghost
! layers are allowed at the physical boundaries - this is consistent
! with our present approach - need to be generalized if only one ghost
! layer is needed
      do iproc=0,numPEs-1
         istart2_all(iproc) = max(imin1-1,min(imax1+1,istart1_all(iproc)-1))
         if(nodesi.ne.1) then
            istart3_all(iproc) = max(imin1-2,min(imax1+2,istart2_all(iproc)-1))
         else
            istart3_all(iproc) = istart2_all(iproc)
         endif

         jstart2_all(iproc) = max(jmin1-1,min(jmax1+1,jstart1_all(iproc)-1))
         if(nodesj.ne.1) then
            jstart3_all(iproc) = max(jmin1-2,min(jmax1+2,jstart2_all(iproc)-1))
         else
            jstart3_all(iproc) = jstart2_all(iproc)
         endif

         kstart2_all(iproc) = max(kmin1-1,min(kmax1+1,kstart1_all(iproc)-1))
         if(nodesk.ne.1) then
            kstart3_all(iproc) = max(kmin1-2,min(kmax1+2,kstart2_all(iproc)-1))
         else
            kstart3_all(iproc) =  kstart2_all(iproc)
         endif

         iend2_all(iproc) = max(imin1-1,min(imax1+1,iend1_all(iproc)+1))
         if(nodesi.ne.1) then
            iend3_all(iproc) = max(imin1-2,min(imax1+2,iend2_all(iproc)+1))
         else
            iend3_all(iproc) = iend2_all(iproc)
         endif

         jend2_all(iproc) = max(jmin1-1,min(jmax1+1,jend1_all(iproc)+1))
         if(nodesj.ne.1) then
            jend3_all(iproc) = max(jmin1-2,min(jmax1+2,jend2_all(iproc)+1))
         else
            jend3_all(iproc) = jend2_all(iproc)
         endif

         kend2_all(iproc) = max(kmin1-1,min(kmax1+1,kend1_all(iproc)+1))
         if(nodesk.ne.1) then
            kend3_all(iproc) = max(kmin1-2,min(kmax1+2,kend2_all(iproc)+1))
         else
            kend3_all(iproc) = kend2_all(iproc)
         endif
      enddo


      do iproc=0,numPEs-1
         ijkstart3_all(iproc) = 1
         ijkend3_all(iproc) =  1 + (iend3_all(iproc) - istart3_all(iproc)) &
           + (jend3_all(iproc)-jstart3_all(iproc))*(iend3_all(iproc)-istart3_all(iproc)+1) &
           + (kend3_all(iproc)-kstart3_all(iproc))*(jend3_all(iproc)-jstart3_all(iproc)+1)* &
             (iend3_all(iproc)-istart3_all(iproc)+1)

      enddo

      do iproc=0,numPEs-1
         ijksize3_all(iproc) = ijkend3_all(iproc) - ijkstart3_all(iproc) + 1
      enddo

      displs(0) = 0
      do iproc=1,numPEs-1
         displs(iproc) = displs(iproc-1)+ijksize3_all(iproc-1)
!       write(*,*) 'displ',displs(iproc),iproc, ijksize3_all(iproc)
      enddo


      ijkstart3 = ijkstart3_all(myPE)
      ijkend3   = ijkend3_all(myPE)
      ijksize3  = ijksize3_all(myPE)

      istart1   = istart1_all(myPE)
      iend1     = iend1_all(myPE)
      jstart1   = jstart1_all(myPE)
      jend1     = jend1_all(myPE)
      kstart1   = kstart1_all(myPE)
      kend1     = kend1_all(myPE)

      istart2   = istart2_all(myPE)
      iend2     = iend2_all(myPE)
      jstart2   = jstart2_all(myPE)
      jend2     = jend2_all(myPE)
      kstart2   = kstart2_all(myPE)
      kend2     = kend2_all(myPE)

      istart3   = istart3_all(myPE)
      iend3     = iend3_all(myPE)
      jstart3   = jstart3_all(myPE)
      jend3     = jend3_all(myPE)
      kstart3   = kstart3_all(myPE)
      kend3     = kend3_all(myPE)


! Setup mapping to take care of cyclic boundary conditions
! ---------------------------------------------------------------->>>
! consider cyclic boundary condition using the imap(:),jmap(:),kmap(:)
! indirection arrays

      allocate( imap( imin4:imax4 ) )
      allocate( jmap( jmin4:jmax4 ) )
      allocate( kmap( kmin4:kmax4 ) )

      allocate( imap_c( imin4:imax4 ) )
      allocate( jmap_c( jmin4:jmax4 ) )
      allocate( kmap_c( kmin4:kmax4 ) )

      do kk=kmin4,kmax4
        kmap(kk) = kk
      enddo

      do jj=jmin4,jmax4
        jmap(jj) = jj
      enddo

      do ii=imin4,imax4
        imap(ii) = ii
      enddo

      if (CYC_ZL) then
         kmap( kmax2 ) = kmin1
         kmap( kmin2 ) = kmax1
         if (kmax3.gt.kmax2) kmap(kmax3) = kmap(kmax2)+1
         if (kmin3.lt.kmin2) kmap(kmin3) = kmap(kmin2)-1
         if (kmax4.gt.kmax3) kmap(kmax4) = kmap(kmax3)+1
         if (kmin4.lt.kmin3) kmap(kmin4) = kmap(kmin3)-1
      endif

      if (CYC_YL) then
         jmap( jmax2 ) = jmin1
         jmap( jmin2 ) = jmax1
         if (jmax3.gt.jmax2) jmap(jmax3) = jmap(jmax2)+1
         if (jmin3.lt.jmin2) jmap(jmin3) = jmap(jmin2)-1
         if (jmax4.gt.jmax3) jmap(jmax4) = jmap(jmax3)+1
         if (jmin4.lt.jmin3) jmap(jmin4) = jmap(jmin3)-1
      endif

      if (CYC_XL) then
         imap( imax2 ) = imin1
         imap( imin2 ) = imax1
         if (imax3.gt.imax2) imap(imax3) = imap(imax2)+1
         if (imin3.lt.imin2) imap(imin3) = imap(imin2)-1
         if (imax4.gt.imax3) imap(imax4) = imap(imax3)+1
         if (imin4.lt.imin3) imap(imin4) = imap(imin3)-1
      endif

      do kk=kmin4,kmax4
        kmap_c(kk) = kk
      enddo

      do jj=jmin4,jmax4
        jmap_c(jj) = jj
      enddo

      do ii=imin4,imax4
         imap_c(ii) = ii
      enddo

      if (CYC_ZL.and.nodesk.eq.1) then
         kmap_c( kmax2 ) = kmin1
         kmap_c( kmin2 ) = kmax1
         if (kmax3.gt.kmax2) kmap_c(kmax3) = kmap_c(kmax2)+1
         if (kmin3.lt.kmin2) kmap_c(kmin3) = kmap_c(kmin2)-1
         if (kmax4.gt.kmax3) kmap_c(kmax4) = kmap_c(kmax3)+1
         if (kmin4.lt.kmin3) kmap_c(kmin4) = kmap_c(kmin3)-1
      endif

      if (CYC_YL.and.nodesj.eq.1) then
         jmap_c( jmax2 ) = jmin1
         jmap_c( jmin2 ) = jmax1
         if (jmax3.gt.jmax2) jmap_c(jmax3) = jmap_c(jmax2)+1
         if (jmin3.lt.jmin2) jmap_c(jmin3) = jmap_c(jmin2)-1
         if (jmax4.gt.jmax3) jmap_c(jmax4) = jmap_c(jmax3)+1
         if (jmin4.lt.jmin3) jmap_c(jmin4) = jmap_c(jmin3)-1
      endif

      if (CYC_XL.and.nodesi.eq.1) then
         imap_c( imax2 ) = imin1
         imap_c( imin2 ) = imax1
         if (imax3.gt.imax2) imap_c(imax3) = imap_c(imax2)+1
         if (imin3.lt.imin2) imap_c(imin3) = imap_c(imin2)-1
         if (imax4.gt.imax3) imap_c(imax4) = imap_c(imax3)+1
         if (imin4.lt.imin3) imap_c(imin4) = imap_c(imin3)-1
      endif
! End setup mapping to take care of cyclic boundary conditions
! ----------------------------------------------------------------<<<


! Defining new set of varaibles to define upper and lower bound of the
! indices to include actual physical boundaries of the problem
      do iproc = 0, numPEs-1
         istart = istart1
         iend = iend1
         jstart = jstart1
         jend = jend1
         kstart = kstart1
         kend = kend1

         if(istart3.eq.imin3) istart = istart2
         if(iend3.eq.imax3) iend = iend2
         if(jstart3.eq.jmin3) jstart = jstart2
         if(jend3.eq.jmax3) jend = jend2
         if(kstart3.eq.kmin3) kstart = kstart2
         if(kend3.eq.kmax3) kend = kend2

         istart_all(iproc) = istart
         iend_all(iproc)   = iend
         jstart_all(iproc) = jstart
         jend_all(iproc)   = jend
         kstart_all(iproc) = kstart
         kend_all(iproc)   = kend
      enddo

      IF(numPEs .GT. 1) THEN
! Calculate any load imbalance.
         IMBALANCE = INT(DBLE(maxval(ijksize3_all)- minval(ijksize3_all)) /    &
            minval(ijksize3_all)*100.0)
! Calculate potential speedup based on Amdahl's Law
         IF(IMBALANCE == 0)THEN
            AMDAHL_SPEEDUP='+Inf'
         ELSE
            AMDAHL_SPEEDUP=''
            WRITE(AMDAHL_SPEEDUP,*)1.0/dble(IMBALANCE)
         ENDIF
! Construct a message for the user telling them the grid partition info.
         WRITE(ERR_MSG,1000)maxval(ijksize3_all), maxloc(ijksize3_all),&
            minval(ijksize3_all), minloc(ijksize3_all),                &
            sum(ijksize3_all)/numPEs, trim(AMDAHL_SPEEDUP)
         CALL FLUSH_ERR_MSG
      ENDIF

! Setup coefficients of FUINIJK
        c0 = 1 - istart3_all(myPE)
        c1 = (iend3_all(myPE)-istart3_all(myPE)+1)
        c2 = (iend3_all(myPE)-istart3_all(myPE)+1)* (jend3_all(myPE)-jstart3_all(myPE)+1)
        c0 =  c0  - c1*jstart3_all(myPE) - c2*kstart3_all(myPE)


! Call to sendrecv_init to set all the communication pattern

      CALL FINL_ERR_MSG
      RETURN

 1000 FORMAT('Parallel load balancing statistics:',2/,13x,'Comp. cells',&
      4X,'Processor',/3X,'maximum   ',I11,4X,I9,/3X,'minimum   ',I11,&
      4X,I9,/3X,'average   ',I11,6X,'-N/A-',2/,3X,'Maximum speedup ',&
      '(Amdahls Law) = ',A)

      END SUBROUTINE GRIDMAP_INIT

      END MODULE GRIDMAP
