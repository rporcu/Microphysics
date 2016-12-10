MODULE READ_RES0_DES_MODULE
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_READ_RESTART                                        !
!  Purpose : Reads either single restart file or multiple restart      !
!  fles (based on bdist_io) flag.                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_RES0_DES(dg_pijk, dg_pijkprv, iglobal_id, particle_state, &
         des_radius, ro_sol, des_usr_var, &
         des_pos_new, des_vel_new, omega_new)

      use des_allocate, only: allocate_dem_mi
      use des_bc, only: dem_mi, dem_bcmi, dem_mi_time
      use discretelement, only: tecplot_findex, des_usr_var_size, vtp_findex, dtsolid
      use error_manager, only: err_msg, flush_err_msg
      use mpi_init_des, only: DES_RESTART_GHOST
      use read_res1_des, only: init_read_res_des, finl_read_res_des, read_par_pos, read_res_des, read_res_parray
      use run, only: run_name, run_type, time

      implicit none

      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: des_radius, ro_sol
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: des_vel_new, des_pos_new, omega_new, des_usr_var
      INTEGER(KIND=1), DIMENSION(:), INTENT(OUT) :: particle_state
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijk, iglobal_id
      INTEGER, DIMENSION(:), INTENT(OUT) :: dg_pijkprv

      INTEGER :: LC1, LC2
      INTEGER :: lDIMN, lNEXT_REC

      DOUBLE PRECISION :: VERSION

      lDIMN = 3

      CALL INIT_READ_RES_DES(trim(RUN_NAME), VERSION, lNEXT_REC)

      CALL READ_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL READ_RES_DES(lNEXT_REC, TECPLOT_FINDEX)
      CALL READ_RES_DES(lNEXT_REC, DTSOLID)

! Position data is read and used to setup pARRAY reads.
      CALL READ_PAR_POS(lNEXT_REC)

      CALL READ_RES_pARRAY(lNEXT_REC, iGLOBAL_ID)

      CALL READ_RES_pARRAY(lNEXT_REC, particle_state)

      DO LC1 = 1, 3
         CALL READ_RES_pARRAY(lNEXT_REC, DES_VEL_NEW(:,LC1))
      ENDDO

      DO LC1 = 1, 3
         CALL READ_RES_pARRAY(lNEXT_REC, OMEGA_NEW(:,LC1))
      ENDDO

      CALL READ_RES_pARRAY(lNEXT_REC, DES_RADIUS)
      CALL READ_RES_pARRAY(lNEXT_REC, RO_SOL)

      IF(VERSION >= 1.1) THEN
         CALL READ_RES_DES(lNEXT_REC, DES_USR_VAR_SIZE)
         DO LC1=1, DES_USR_VAR_SIZE
            CALL READ_RES_pARRAY(lNEXT_REC, DES_USR_VAR(:,LC1))
         ENDDO
      ENDIF

      CALL DES_RESTART_GHOST(dg_pijk, dg_pijkprv, des_pos_new)

! RES2 does not need the collision of BC information.
      IF(RUN_TYPE == 'RESTART_2') RETURN

! Collision/neighbor data is read and used to setup cARRAY reads.
!      CALL READ_PAR_COL(lNEXT_REC)
!
!      DO LC1=1, lDIMN
!         CALL READ_RES_cARRAY(lNEXT_REC, PFT_NEIGHBOR(LC1,:))
!      ENDDO

! Save the number of BCMI's read from input file, then read the
! value from the restart file.
      CALL READ_RES_DES(lNEXT_REC, DEM_BCMI)

! Allocation of MIs is done here to ignore changes to the mfix.dat
! file during RES1.
      IF(DEM_BCMI > 0) CALL ALLOCATE_DEM_MI

! Only save the number of mass inflows for RESTART_1. This allows
! for mass inflows to be added/removed with RESTART_2.
! Todo: Prune entering/exiting flagged particles for RESTART_2.
      DO LC1=1, DEM_BCMI
         CALL READ_RES_DES(lNEXT_REC, DEM_MI_TIME(LC1))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%VACANCY)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%OCCUPANTS)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%WINDOW)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%OFFSET)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%L)

         LC2 = DEM_MI(LC1)%OCCUPANTS

         allocate(DEM_MI(LC1)%W(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%W(:))
         allocate(DEM_MI(LC1)%H(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%H(:))
         allocate(DEM_MI(LC1)%P(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%P(:))
         allocate(DEM_MI(LC1)%Q(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%Q(:))
      ENDDO

      CALL FINL_READ_RES_DES


      WRITE(ERR_MSG,"('DES restart file read at Time = ',g12.5)") TIME
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      RETURN
      END SUBROUTINE READ_RES0_DES
END MODULE READ_RES0_DES_MODULE
