!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_OUTFLOW(BCV, I1, I2, J1, J2, K1, K2)               C
!  Author: M. Syamlal                                 Date: 21-JAN-92  C
!                                                                      C
!  Purpose: Set specified pressure outflow bc for a specified range of C
!           cells                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CG_SET_OUTFLOW

      USE bc
      USE compar        !//d
      USE cutcell
      USE eos, ONLY: EOSG
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE mflux
      USE param
      USE param1
      USE physprop
      USE quadric
      USE run

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      indices
      INTEGER          I, J, K
!
!                      Local index for boundary cell
      INTEGER          IJK
!
!                      Boundary condition number
      INTEGER          BCV
!
!                      Locall index for a fluid cell near the boundary cell
      INTEGER          LFLUID

      INTEGER          IJKW,IJKWW,IJKS,IJKSS,IJKB

      CHARACTER(LEN=16) :: BCT1,BCT2,BCT3,BCT4

      LOGICAL :: TEST1,TEST2
!-----------------------------------------------
!
!      print*,'top of cg_set_outflow'
      DO IJK = IJKSTART3, IJKEND3
      IF(INTERIOR_CELL_AT(IJK)) THEN


         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

!//SP Check if current i,j,k resides on this PE
         IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
!        IJK = FUNIJK(I,J,K)
!
! Fluid cell at West
!
!      print*,'west'
          BCV = BC_U_ID(IJK)


          IF(BCV>0)THEN

            IF (BC_TYPE(BCV) == 'CG_PO') THEN

               IF (FLUID_AT(IM_OF(IJK))) THEN
                  LFLUID = IM_OF(IJK)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(mw_avg,P_G&
                        (IJK),295.15d0)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF


            ENDIF
         ENDIF

         IJKW = WEST_OF(IJK)
         IJKWW = WEST_OF(IJKW)

         BCT1=''
         BCT2=''
         BCT3=''
         BCV = BC_ID(IJK)
         IF(BCV>0) BCT1 = BC_TYPE(BCV)
         BCV = BC_U_ID(IJK)
         IF(BCV>0) BCT2 = BC_TYPE(BCV)
         BCV = BC_U_ID(IJKW)
         IF(BCV>0) BCT3 = BC_TYPE(BCV)
         BCV = BC_ID(IJKW)
         IF(BCV>0) BCT4 = BC_TYPE(BCV)

         TEST1= ((BCT1 == 'CG_PO').AND.&
            (BCT2 /= 'CG_PO').AND.(BCT3 == 'CG_PO'))

         TEST2 = (BCT2 == 'CG_PO').AND.(BCT4 /= 'CG_PO')



         IF(TEST1) THEN
            BCV = BC_ID(IJK)

            RO_G(IJK) = RO_G(IJKW)
            EP_G(IJK) = EP_G(IJKW)
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)

         ENDIF




!
! Fluid cell at East
!      print*,'east'

          BCV = BC_U_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE(BCV) == 'CG_PO') THEN


               IF (FLUID_AT(IP_OF(IJK))) THEN
                  LFLUID = IP_OF(IJK)
                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(mw_avg,P_G&
                        (IJK),295.15d0)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF

            ENDIF
         ENDIF

!
! Fluid cell at South
!
!      print*,'south'
          BCV = BC_V_ID(IJK)
!      print*,ijk,I,J,K,bcv,BC_TYPE(BCV)
          IF(BCV>0)THEN

            IF (BC_TYPE(BCV) == 'CG_PO') THEN


               IF (FLUID_AT(JM_OF(IJK))) THEN
                  LFLUID = JM_OF(IJK)

                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(mw_avg,P_G&
                        (IJK),295.15d0)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE

                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF

            ENDIF
         ENDIF

         IJKS = SOUTH_OF(IJK)
         IJKSS = SOUTH_OF(IJKS)
!         print*,'ijk,ijks=',ijk,ijks,JM_OF(IJK)
         BCT1=''
         BCT2=''
         BCT3=''
         BCV = BC_ID(IJK)
!         print*,'bcv=',bcv
         IF(BCV>0) BCT1 = BC_TYPE(BCV)
         BCV = BC_V_ID(IJK)
!         print*,'bcv=',bcv
         IF(BCV>0) BCT2 = BC_TYPE(BCV)
         BCV = BC_V_ID(IJKS)
!         print*,'bcv=',bcv
         IF(BCV>0) BCT3 = BC_TYPE(BCV)

!         IF((BCT1 == 'CG_PO').AND.(BCT2 /= 'CG_PO').AND.(BCT3 == 'CG_PO')) THEN



         TEST1 = (BCT1 == 'CG_PO').AND.&
            (BCT2 /= 'CG_PO').AND.(BCT3 == 'CG_PO')
         TEST2 = (FLAG(IJK) == 11).AND.(BCT2 == 'CG_PO')

         TEST2 = (FLAG(IJK) == 11)

         IF(TEST1) THEN

            BCV = BC_ID(IJK)

            RO_G(IJK) = RO_G(IJKS)
            EP_G(IJK) = EP_G(IJKS)
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ENDIF



!
! Fluid cell at North
!

!      print*,'north'
          BCV = BC_V_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE(BCV) == 'CG_PO') THEN


               IF (FLUID_AT(JP_OF(IJK))) THEN
                  LFLUID = JP_OF(IJK)

                     IF (RO_G0 == UNDEFINED) RO_G(IJK) = EOSG(mw_avg,P_G&
                        (IJK),29.15d0)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF

            ENDIF
         ENDIF

!
! Fluid cell at Bottom
!
!      print*,'bottom'

          BCV = BC_W_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE(BCV) == 'CG_PO') THEN

               IF (FLUID_AT(KM_OF(IJK))) THEN
                  LFLUID = KM_OF(IJK)
                  IF (RO_G0 == UNDEFINED) RO_G(IJK) = &
                     EOSG(mw_avg, P_G(IJK),29.15d0)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF


            ENDIF
         ENDIF

         IJKB = BOTTOM_OF(IJK)
         BCT1=''
         BCT2=''
         BCT3=''
         BCV = BC_ID(IJK)
         IF(BCV>0) BCT1 = BC_TYPE(BCV)
         BCV = BC_W_ID(IJK)
         IF(BCV>0) BCT2 = BC_TYPE(BCV)
         BCV = BC_W_ID(IJKB)
         IF(BCV>0) BCT3 = BC_TYPE(BCV)

         IF((BCT1 == 'CG_PO').AND.(BCT2 /= 'CG_PO').AND.(BCT3 == 'CG_PO')) THEN
!            print*,'IJK,IJKB=',IJK,IJKB
!            read(*,*)
            BCV = BC_ID(IJK)
            RO_G(IJK) = RO_G(IJKB)
            EP_G(IJK) = EP_G(IJKB)
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ENDIF


!
! Fluid cell at Top
!
!      print*,'top'
          BCV = BC_W_ID(IJK)

          IF(BCV>0)THEN

            IF (BC_TYPE(BCV) == 'CG_PO') THEN


               IF (FLUID_AT(KP_OF(IJK))) THEN
                  LFLUID = KP_OF(IJK)
                  IF (RO_G0 == UNDEFINED) RO_G(IJK) = &
                     EOSG(mw_avg,P_G(IJK),295.15d0)
                  IF (BC_EP_G(BCV) == UNDEFINED) EP_G(IJK) = ONE
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF

            ENDIF
         ENDIF

      ENDIF
      END DO

!      print*,'bottom of cg_set_outflow'
      RETURN
      END SUBROUTINE CG_SET_OUTFLOW
