module comp_mean_fields0_module

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS0(ep_g,ro_g,rop_g,pijk,particle_phase,pmass,pvol,des_pos_new,des_vel_new)

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar, only: iend1, jend1, kend1
      USE compar, only: imap_c, jmap_c, kmap_c
      USE compar, only: istart2, jstart2, kstart2
      USE compar, only: mype, pe_io
      USE discretelement, only: des_rop_s, des_rops_node, xe, yn, zt, interp_scheme
      USE discretelement, only: pic, des_vel_node, dimn, pinc
      USE calc_epg_des_module, only: calc_epg_des
      USE interpolation, only: set_interpolation_scheme
      USE functions, only: fluid_at
      USE geometry, only: vol_surr, vol
      USE interpolation, only: set_interpolation_stencil
      USE mpi_node_des, only: des_addnodevalues_mean_fields
      USE param1, only: zero
      USE particle_filter, only: DES_REPORT_MASS_INTERP
      USE constant, only: mmax, ro_s0

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ep_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(IN   ) :: ro_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      DOUBLE PRECISION, INTENT(INOUT) :: rop_g&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pmass, pvol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new
      INTEGER, DIMENSION(:,:), INTENT(IN) :: pijk
      INTEGER, DIMENSION(:), INTENT(IN) :: particle_phase

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, &
                 II, JJ, KK
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3):: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to
! set_interpolation_stencil
      INTEGER :: ONEW
! constant whose value depends on dimension of system
! avg_factor=0.250 (in 3D) or =0.50 (in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
! index of solid phase that particle NP belongs to
      INTEGER :: llI, llJ, llK, M
! particle number index, used for looping
      INTEGER :: NP, NINDX

      DOUBLE PRECISION :: MASS_SOL1, MASS_SOL2
! sum of mass_sol1 and mass_sol2 across all processors
      DOUBLE PRECISION :: MASS_SOL1_ALL, MASS_SOL2_ALL

      DOUBLE PRECISION :: TEMP1, TEMP2

      DOUBLE PRECISION, DIMENSION(3) :: DES_VEL_DENSITY
      DOUBLE PRECISION :: DES_ROP_DENSITY

      INTEGER :: COUNT_NODES_OUTSIDE, &
                 COUNT_NODES_INSIDE_MAX
!Handan Liu added on Jan 17 2013
          DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp
          DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft
!-----------------------------------------------

! initializing
      MASS_SOL1 = ZERO
      MASS_SOL2 = ZERO
      MASS_SOL1_ALL = ZERO
      MASS_SOL2_ALL = ZERO
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = 0.25D0

! cartesian_grid related quantities
      COUNT_NODES_INSIDE_MAX =8


! Initialize entire arrays to zero
      DES_VEL_NODE = ZERO
      DES_ROPS_NODE = ZERO
      DES_ROP_S = zero


! sets several quantities including interp_scheme, scheme, and
! order and allocates arrays necessary for interpolation
      CALL SET_INTERPOLATION_SCHEME(2)

      DO llK = kstart3, kend3
         DO llJ = jstart3, jend3
            DO llI = istart3, iend3

! Cycle this cell if not in the fluid domain or if it contains no
! particle/parcel
               IF(.NOT.fluid_at(lli,llj,llk)) CYCLE
               IF( PINC(lli,llj,llk) == 0) CYCLE

               PCELL(1) = lli-1
               PCELL(2) = llj-1
               PCELL(3) = llk-1

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
               CALL SET_INTERPOLATION_STENCIL(PCELL, IW, IE, JS, JN, KB, KTP,&
                  INTERP_SCHEME, DIMN, ORDERNEW=ONEW)

               COUNT_NODES_OUTSIDE = 0
! Computing/setting the geometric stencil
               DO K=1, ONEW
                  DO J=1, ONEW
                     DO I=1, ONEW

                        GST_TMP(I,J,K,1) = XE(IW + I-1)
                        GST_TMP(I,J,K,2) = YN(JS + J-1)
                        GST_TMP(I,J,K,3) = ZT(KB + K-1)

                     ENDDO
                  ENDDO
               ENDDO


! Calculate des_rops_node so des_rop_s, and in turn, ep_g can be updated
!----------------------------------------------------------------->>>

! looping through particles in the cell
               DO NINDX=1, PINC(lli,llj,llk)
                  NP = PIC(lli,llj,llk)%P(NINDX)

                  call DRAG_WEIGHTFACTOR(gst_tmp,des_pos_new(np,:),weight_ft)

                  M = particle_phase(NP)

                  MASS_SOL1 = MASS_SOL1 + PMASS(NP)

                  TEMP2 = RO_S0(M)*PVOL(NP)

                  DO K = 1, ONEW
                     DO J = 1, ONEW
                        DO I = 1, ONEW
! shift loop index to new variables for manipulation
                           II = IW + I-1
                           JJ = JS + J-1
                           KK = KB + K-1

                           TEMP1 = WEIGHT_FT(I,J,K)*TEMP2

                           DES_ROPS_NODE(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK),M) = &
                              DES_ROPS_NODE(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK),M) + TEMP1

                           DES_VEL_NODE(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK),:,M) = &
                              DES_VEL_NODE(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK),:,M) + &
                              TEMP1*DES_VEL_NEW(NP,:)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO   ! end do (nindx=1,pinc(ijk))
!-----------------------------------------------------------------<<<

            ENDDO
         ENDDO
      ENDDO


! At the interface des_rops_node has to be added since particles
! across the processors will contribute to the same scalar node.
! sendrecv will be called and the node values will be added
! at the junction. des_rops_node is altered by the routine when
! periodic boundaries are invoked
      CALL DES_ADDNODEVALUES_MEAN_FIELDS


! Now go from node to scalar center. Same convention as sketched
! earlier
!----------------------------------------------------------------->>>
! Explanation by RG: 08/17/2012
! the approach used here is to make it general enough for cutcells to be
! included as well. The new changes do not alter earlier calculations
! but make the technique general as to include cartesian grid (cut-cell)
! simulations.
! Previously, the volume of the node (by array des_vol_node) was used to
! first scale the nodal values. Subsequently, these nodal values were
! equally added to compute the cell centered values for the scalar cell.

! Consider an internal node next to an edge node (a node adjacent to a
! boundary). In 2D, the volume of an edge node will be half that of an
! internal node. And an edge node will contribute double compared to
! an internal node to the value of the scalar cell they share. These
! calculations were previously accomplished via the variable volume of
! node.  Now this is accomplished by the ratio vol(ijk2)/vol_sur, where
! vol(ijk2) is the volume of the scalar cell in consideration and
! vol_sur is the sum of all the scalar cell volumes that have this node
! as the common node.
!---------------------------------------------------------------------//
      DO K = KSTART2, KEND1
      DO J = JSTART2, JEND1
      DO I = ISTART2, IEND1
         if (vol_surr(i,j,k).eq.ZERO) CYCLE

! looping over stencil points (NODE VALUES)
         DO M = 1, MMAX

            DES_ROP_DENSITY = DES_ROPS_NODE(i,j,k, M)/vol_surr(i,j,k)
            DES_VEL_DENSITY(:) = DES_VEL_NODE(i,j,k, :, M)/vol_surr(i,j,k)

            DO KK = K, K+1
            DO JJ = J, J+1
            DO II = I, I+1

               IF(fluid_at(II,JJ,KK)) THEN
! Since the data in the ghost cells is spurious anyway and overwritten during
! subsequent send receives, do not compute any value here as this will
! mess up the total mass value that is computed below to ensure mass conservation
! between Lagrangian and continuum representations
                  DES_ROP_S(IMAP_C(ii), JMAP_C(jj), KMAP_C(kk), M) = &
                     DES_ROP_S(IMAP_C(ii), JMAP_C(jj), KMAP_C(kk), M) + &
                     DES_ROP_DENSITY*VOL
               ENDIF
            ENDDO  ! end do (ii=i1,i2)
            ENDDO  ! end do (jj=j1,j2)
            ENDDO  ! end do (kk=k1,k2)
         ENDDO
      ENDDO   ! end do (i=istart2,iend1)
      ENDDO   ! end do (j=jstart2,jend1)
      ENDDO   ! end do (k=kstart2,kend1)

!-----------------------------------------------------------------<<<


      DO llK = kstart3, kend3
      DO llJ = jstart3, jend3
      DO llI = istart3, iend3
         IF(.NOT.fluid_at(lli,llj,llk)) CYCLE

         DO M = 1, MMAX
            IF(DES_ROP_S(lli,llj,llk, M).GT.ZERO) THEN

! Divide by scalar cell volume to obtain the bulk density
               DES_ROP_S(lli,llj,llk, M) = DES_ROP_S(lli,llj,llk, M)/VOL

            ENDIF
         ENDDO   ! end loop over M=1,MMAX
      ENDDO
      ENDDO
      ENDDO

      CALL CALC_EPG_DES(ep_g,ro_g,rop_g,des_pos_new)

! turn on the below statements to check if the mass is conserved
! between discrete and continuum representations. Should be turned to
! false for any production runs.
      IF(DES_REPORT_MASS_INTERP) THEN

      DO llK = kstart3, kend3
        DO llJ = jstart3, jend3
         DO llI = istart3, iend3

! It is important to check fluid_at
            IF(.NOT.fluid_at(lli,llj,llk)) CYCLE

            MASS_SOL2 = MASS_SOL2 + sum(DES_ROP_S(lli,llj,llk,1:MMAX))*VOL
         ENDDO
         ENDDO
         ENDDO


         mass_sol1_all = mass_sol1
         ! CALL GLOBAL_SUM(MASS_SOL1, MASS_SOL1_ALL)
         mass_sol2_all = mass_sol2
         ! CALL GLOBAL_SUM(MASS_SOL2, MASS_SOL2_ALL)
         if(myPE.eq.pe_IO) THEN
            WRITE(*,'(/,5x,A,4(2x,g17.8),/)') &
                 'SOLIDS MASS DISCRETE AND CONTINUUM =  ', &
                 MASS_SOL1_ALL, MASS_SOL2_ALL
         ENDIF
      ENDIF
      END SUBROUTINE COMP_MEAN_FIELDS0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!  Subroutine: DRAG_WEIGHTFACTOR                                        C
!  Purpose: DES - Calculate the fluid velocity interpolated at the      C
!           particle's location and weights. Replace 'interpolator'     C
!                       interface for OpenMP implementation.            C
!                                                                       C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE DRAG_WEIGHTFACTOR(GSTEN,DESPOS,WEIGHTFACTOR)

        IMPLICIT NONE

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
        DOUBLE PRECISION, DIMENSION(2,2,2,3), INTENT(IN):: GSTEN
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: DESPOS
        DOUBLE PRECISION, DIMENSION(2,2,2), INTENT(OUT) :: WEIGHTFACTOR
        INTEGER :: II, JJ, KK

        DOUBLE PRECISION, DIMENSION(2) :: XXVAL, YYVAL, ZZVAL
        DOUBLE PRECISION :: DXX, DYY, DZZ
        DOUBLE PRECISION, DIMENSION(3) :: ZETAA

        DXX = GSTEN(2,1,1,1) - GSTEN(1,1,1,1)
        DYY = GSTEN(1,2,1,2) - GSTEN(1,1,1,2)

        ZETAA(1:2) = DESPOS(1:2) - GSTEN(1,1,1,1:2)

        ZETAA(1) = ZETAA(1)/DXX
        ZETAA(2) = ZETAA(2)/DYY

        XXVAL(1)=1-ZETAA(1)
        YYVAL(1)=1-ZETAA(2)
        XXVAL(2)=ZETAA(1)
        YYVAL(2)=ZETAA(2)

        DZZ = GSTEN(1,1,2,3) - GSTEN(1,1,1,3)
        ZETAA(3) = DESPOS(3) - GSTEN(1,1,1,3)
        ZETAA(3) = ZETAA(3)/DZZ
        ZZVAL(1)=1-ZETAA(3)
        ZZVAL(2)=ZETAA(3)
        DO KK=1,2
           DO JJ=1,2
              DO II=1,2
                 WEIGHTFACTOR(II,JJ,KK) = XXVAL(II)*YYVAL(JJ)*ZZVAL(KK)
              ENDDO
           ENDDO
        ENDDO

      END SUBROUTINE DRAG_WEIGHTFACTOR

end module comp_mean_fields0_module
