MODULE WRITE_RES0_DES_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
!
!  module name: des_write_restart
!  purpose: writing des data for restart
!
!  Author : Pradeep G
!  Purpose : Reads either single restart file or multiple restart files
!            (based on bdist_io) flag
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^c
      SUBROUTINE WRITE_RES0_DES(iglobal_id, particle_state, &
         des_radius, ro_sol, des_usr_var, &
         des_pos_new, des_vel_new, omega_new)

      use discretelement, only: des_usr_var_size, dtsolid, tecplot_findex, vtp_findex
      use run, only: run_name
      use des_bc, only: dem_mi, dem_bcmi, dem_mi_time

      use write_res1_des, only: init_write_res_des, write_res_des, write_res_parray, finl_write_res_des

      implicit none

      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: des_radius, ro_sol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: des_vel_new, des_pos_new, omega_new, des_usr_var
      INTEGER(KIND=1), DIMENSION(:), INTENT(IN) :: particle_state
      INTEGER, DIMENSION(:), INTENT(IN) :: iglobal_id

!-----------------------------------------------
! local variables
!-----------------------------------------------
      INTEGER :: LC1
      INTEGER :: lNEXT_REC
      INTEGER :: lDIMN

      DOUBLE PRECISION :: VERSION

! Set the version of the DES RES file.
      VERSION = 1.1

! Set the output dimension.
      lDIMN = 3

      CALL INIT_WRITE_RES_DES(trim(RUN_NAME), VERSION, lNEXT_REC)

      CALL WRITE_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, TECPLOT_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, DTSOLID)

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_POS_NEW(:,LC1))
      ENDDO

      CALL WRITE_RES_pARRAY(lNEXT_REC, iGLOBAL_ID)

      CALL WRITE_RES_pARRAY(lNEXT_REC, particle_state)

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_VEL_NEW(:,LC1))
      ENDDO

      DO LC1 = 1, 3
         CALL WRITE_RES_pARRAY(lNEXT_REC, OMEGA_NEW(:,LC1))
      ENDDO

      CALL WRITE_RES_pARRAY(lNEXT_REC, DES_RADIUS)
      CALL WRITE_RES_pARRAY(lNEXT_REC, RO_SOL)

! DES User defined variable :: added for VERSION>= 1.1
      CALL WRITE_RES_DES(lNEXT_REC, DES_USR_VAR_SIZE)
      DO LC1=1,DES_USR_VAR_SIZE
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_USR_VAR(:,LC1))
      ENDDO

!      CALL WRITE_RES_cARRAY(lNEXT_REC, NEIGHBORS(:), pLOC2GLB=.TRUE.)
!      CALL WRITE_RES_cARRAY(lNEXT_REC, NEIGHBOR_INDEX(:), pLOC2GLB=.TRUE.)
!
!      DO LC1=1, lDIMN
!         CALL WRITE_RES_cARRAY(lNEXT_REC,PFT_NEIGHBOR(LC1,:))
!      ENDDO

      CALL WRITE_RES_DES(lNEXT_REC, DEM_BCMI)
      DO LC1=1, DEM_BCMI
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI_TIME(LC1))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%VACANCY)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OCCUPANTS)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%WINDOW)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OFFSET)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%L)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%W(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%H(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%P(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%Q(:))
      ENDDO

      CALL FINL_WRITE_RES_DES

      RETURN
      END SUBROUTINE WRITE_RES0_DES
END MODULE WRITE_RES0_DES_MODULE
