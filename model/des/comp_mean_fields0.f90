module comp_mean_fields0_module

   contains
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE COMP_MEAN_FIELDS0(ep_g,ro_g,rop_g,particle_phase,pmass,pvol, &
         des_pos_new,des_vel_new,des_radius,des_usr_var,vol_surr,iglobal_id,flag,pinc)

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
      USE compar, only: iend1, jend1, kend1
      USE compar, only: imap_c, jmap_c, kmap_c
      USE compar, only: istart2, jstart2, kstart2
      USE compar, only: mype, pe_io
      USE constant, only: mmax
      USE discretelement, only: des_rops_node, xe, yn, zt, interp_scheme
      USE discretelement, only: pic, dimn
      USE interpolation, only: set_interpolation_scheme
      USE geometry, only: vol
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
      INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: FLAG
      DOUBLE PRECISION, INTENT(in   ) :: vol_surr&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      INTEGER,          INTENT(inout) :: pinc&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)

      double precision, intent(in) :: des_radius(:)
      double precision, intent(in) :: pmass(:), pvol(:)
      double precision, intent(in) :: des_vel_new(:,:)
      double precision, intent(in) :: des_pos_new(:,:)
      double precision, intent(in) :: des_usr_var(:,:)
      integer, intent(in) :: particle_phase(:)
      integer, intent(in) :: iglobal_id(:)

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
! index of solid phase that particle NP belongs to
      INTEGER :: llI, llJ, llK
! particle number index, used for looping
      INTEGER :: NP, NINDX

      DOUBLE PRECISION :: TEMP1, TEMP2

      DOUBLE PRECISION :: gst_tmp(2,2,2,3)
      DOUBLE PRECISION :: weight_ft(2,2,2)

      double precision :: OoVOL

      double precision :: PVOL_NODE&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
      double precision :: EPs&
         (istart3:iend3, jstart3:jend3, kstart3:kend3)
!-----------------------------------------------

! Initialize entire arrays to zero
      PVOL_NODE = 0.0d0
      EPs = 0.0d0

! sets several quantities including interp_scheme, scheme, and
! order and allocates arrays necessary for interpolation
      CALL SET_INTERPOLATION_SCHEME(2)

      DO llK = kstart3, kend3
         DO llJ = jstart3, jend3
            DO llI = istart3, iend3

! Cycle this cell if not in the fluid domain or if it contains no
! particle/parcel
               IF(.NOT.1.eq.flag(lli,llj,llk,1)) CYCLE
               IF( PINC(lli,llj,llk) == 0) CYCLE

               PCELL(1) = lli-1
               PCELL(2) = llj-1
               PCELL(3) = llk-1

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
               CALL SET_INTERPOLATION_STENCIL(PCELL, IW, IE, JS, JN, KB, KTP,&
                  INTERP_SCHEME, DIMN, ORDERNEW=ONEW)

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

!----------------------------------------------------------------->>>

! looping through particles in the cell
               DO NINDX=1, PINC(lli,llj,llk)
                  NP = PIC(lli,llj,llk)%P(NINDX)

                  call DRAG_WEIGHTFACTOR(gst_tmp,des_pos_new(np,:),weight_ft)

                  DO K = 1, ONEW
                     DO J = 1, ONEW
                        DO I = 1, ONEW
! shift loop index to new variables for manipulation
                           II = IW + I-1
                           JJ = JS + J-1
                           KK = KB + K-1

                           TEMP1 = WEIGHT_FT(I,J,K)*PVOL(NP)

                           PVOL_NODE(II,JJ,KK) = PVOL_NODE(II,JJ,KK) + TEMP1

                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO   ! end do (nindx=1,pinc(ijk))
!-----------------------------------------------------------------<<<

            ENDDO
         ENDDO
      ENDDO



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

               TEMP1 = VOL*(PVOL_NODE(i,j,k)/vol_surr(i,j,k))

               DO KK = K, K+1
                  DO JJ = J, J+1
                     DO II = I, I+1

                        IF(1.eq.flag(II,JJ,KK,1)) THEN
                           EPs(ii,jj,kk) = EPs(ii,jj,kk) + TEMP1
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

!-----------------------------------------------------------------<<<

      OoVOL = 1.0d0/VOL
      DO llK = kstart3, kend3
         DO llJ = jstart3, jend3
            DO llI = istart3, iend3
               IF(.NOT.1.eq.flag(lli,llj,llk,1)) CYCLE

! ! Divide by scalar cell volume to obtain the bulk density
               ep_g(lli,llj,llk) = 1.0d0 - EPs(lli,llj,llk)*OoVOL
            ENDDO
         ENDDO
      ENDDO



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
