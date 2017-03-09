!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: bc.inc                                                 C
!  Purpose: Common block containing boundary conditions data           C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: None                                C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE bc

      use amrex_fort_module, only : c_real => amrex_real
      use iso_c_binding , only: c_int

      USE param, only: dimension_bc, dim_m, dim_n_g, dim_n_s, dimension_ic
      USE param1, only: zero, small_number, one, undefined, half

!                      x coordinate of the west face of a region where
!                      boundary conditions are specified
      real(c_real) BC_X_w (DIMENSION_BC)
!
!                      x coordinate of the east face of a region where
!                      boundary conditions are specified
      real(c_real) BC_X_e (DIMENSION_BC)
!
!                      y coordinate of the south face of a region where
!                      boundary conditions are specified
      real(c_real) BC_Y_s (DIMENSION_BC)
!
!                      y coordinate of the north face of a region where
!                      boundary conditions are specified
      real(c_real) BC_Y_n (DIMENSION_BC)
!
!                      z coordinate of the bottom face of a region where
!                      boundary conditions are specified
      real(c_real) BC_Z_b (DIMENSION_BC)
!
!                      z coordinate of the top face of a region where
!                      boundary conditions are specified
      real(c_real) BC_Z_t (DIMENSION_BC)
!
!                      i index of the west face of a region where
!                      boundary conditions are specified
      INTEGER          BC_I_w (DIMENSION_BC)
!
!                      i index of the east face of a region where
!                      boundary conditions are specified
      INTEGER          BC_I_e (DIMENSION_BC)
!
!                      j index of the south face of a region where
!                      boundary conditions are specified
      INTEGER          BC_J_s (DIMENSION_BC)
!
!                      j index of the north face of a region where
!                      boundary conditions are specified
      INTEGER          BC_J_n (DIMENSION_BC)
!
!                      k index of the bottom face of a region where
!                      boundary conditions are specified
      INTEGER          BC_K_b (DIMENSION_BC)
!
!                      k index of the top face of a region where
!                      boundary conditions are specified
      INTEGER          BC_K_t (DIMENSION_BC)
!
!                      Void fraction in a specified boundary region
      real(c_real) BC_EP_g (DIMENSION_BC)
!
!                      Gas pressure in a specified boundary region
      real(c_real) BC_P_g (DIMENSION_BC)
!
!                      Microscopic density of gas in a specified
!                      boundary region
      real(c_real) BC_RO_g (DIMENSION_BC)
!
!                      Macroscopic density of gas in a specified
!                      boundary region
      real(c_real) BC_ROP_g (DIMENSION_BC)
!
!                      Macroscopic density of solids phases in a
!                      specified boundary region
      real(c_real) BC_ROP_s (DIMENSION_BC, DIM_M)
      real(c_real) BC_EP_s (DIMENSION_BC, DIM_M)
!
!                      Gas phase temperature in a specified boundary
!                      region
      real(c_real) BC_T_g (DIMENSION_BC)
!
!                      Solids phase temperature in a specified
!                      boundary region
      real(c_real) BC_T_s (DIMENSION_BC, DIM_M)
!
!
!                      x-component of gas velocity in a specified
!                      boundary region
      real(c_real) BC_U_g (DIMENSION_BC)
!
!                      x-component of solids phase velocity in a
!                      specified boundary region
      real(c_real) BC_U_s (DIMENSION_BC, DIM_M)
!
!                      y-component of gas velocity in a specified
!                      boundary region
      real(c_real) BC_V_g (DIMENSION_BC)
!
!                      y-component of solids phase velocity in a
!                      specified boundary region
      real(c_real) BC_V_s (DIMENSION_BC, DIM_M)
!
!                      z-component of gas velocity in a specified
!                      boundary region
      real(c_real) BC_W_g (DIMENSION_BC)
!
!                      z-component of solids phase velocity in a
!                      specified boundary region
      real(c_real) BC_W_s (DIMENSION_BC, DIM_M)
!
! JFD: For cut cells, define the magnitude of velocity that will be enforced
!      perpendicular to the cut face, for CG_MI boundary condition
!
!                      magnitude of gas velocity in a specified
!                      boundary region
      real(c_real) BC_VELMAG_g (DIMENSION_BC)
!
!                      magnitude of solids phase velocity in a
!                      specified boundary region
      real(c_real) BC_VELMAG_s (DIMENSION_BC, DIM_M)

!
!                      Type of boundary: MASS_INFLOW, MASS_OUTFLOW,
!                      P_INFLOW, P_OUTFLOW, FREE_SLIP_WALL, NO_SLIP_WALL
      CHARACTER(LEN=16)     BC_TYPE (DIMENSION_BC)
      INTEGER ::     BC_TYPE_ENUM (DIMENSION_BC)

      ENUM, BIND(C)
         ENUMERATOR :: CG_NSW, CG_FSW, CG_PSW, CG_MI, NONE
         ENUMERATOR :: NO_SLIP_WALL, FREE_SLIP_WALL, PAR_SLIP_WALL, NSW, FSW, PSW
         ENUMERATOR :: P_OUTFLOW, MASS_OUTFLOW, OUTFLOW
         ENUMERATOR :: P_INFLOW, MASS_INFLOW
         ENUMERATOR :: CG_PO, CG_MO
         ENUMERATOR :: BLANK
         ENUMERATOR :: DUMMY
      END ENUM

!                      FLAG to specify if this PO BC applies to solid phase
!                      in discrete implementation or not. For example, setting
!                      Pressure outflow only for gas-phase.

      LOGICAL          BC_PO_APPLY_TO_DES (DIMENSION_BC)

!                      Gas volumetric flow rate through the boundary
      real(c_real) BC_VOLFLOW_g (DIMENSION_BC)
!
!                      Solids volumetric flow rate through the boundary
      real(c_real) BC_VOLFLOW_s (DIMENSION_BC, DIM_M)
!
!                      Gas mass flow rate through the boundary
      real(c_real) BC_MASSFLOW_g (DIMENSION_BC)
!
!                      Solids mass flow rate through the boundary
      real(c_real) BC_MASSFLOW_s (DIMENSION_BC, DIM_M)
!
!                      Logical variable to determine whether a bc is defined
       LOGICAL         BC_DEFINED (DIMENSION_BC)
!
!                      Character variable with values W, E, S, N, B, and T
!                      to determine the flow plane of a flow cell
       CHARACTER       bc_plane (DIMENSION_BC)
!
!                      The interval at the beginning when normal vel. is equal to
!                      BC_Jet_g0
      real(c_real) BC_DT_0 (DIMENSION_BC)
!
!                      Stored value of normal velocity
      real(c_real) BC_Jet_g (DIMENSION_BC)
!
!                      Value of normal vel. during the initial interval BC_DT_0
      real(c_real) BC_Jet_g0 (DIMENSION_BC)
!
!                      The interval when normal vel. is equal to BC_Jet_gh
      real(c_real) BC_DT_h (DIMENSION_BC)
!
!                      Value of normal vel. during the initial interval BC_DT_h
      real(c_real) BC_Jet_gh (DIMENSION_BC)
!
!                      The interval when normal vel. is equal to BC_Jet_gl
      real(c_real) BC_DT_l (DIMENSION_BC)
!
!                      Value of normal vel. during the initial interval BC_DT_l
      real(c_real) BC_Jet_gl (DIMENSION_BC)
!
!                      Time to update a transient boundary condition
      real(c_real) BC_TIME (DIMENSION_BC)
!
!                      Area of boundary surfaces
      real(c_real) BC_AREA (DIMENSION_BC)
!
!                      Volume of boundary cells
      real(c_real) BC_VOL (DIMENSION_BC)
!
!                      Gas species mass fractions in a boundary region
      real(c_real) BC_X_g (DIMENSION_BC, DIM_N_g)
!
!                      Solids species mass fractions in a boundary region
      real(c_real) BC_X_s (DIMENSION_BC, DIM_M, DIM_N_s)
!
!                      Accumulated or average mass outflow rate of gas
      real(c_real) BC_MOUT_g(DIMENSION_BC)
!
!                      Accumulated or average mass outflow rate of solids
      real(c_real) BC_MOUT_s(DIMENSION_BC, DIM_M)
!
!                      Accumulated or average volumetric outflow rate of gas
      real(c_real) BC_VOUT_g(DIMENSION_BC)
!
!                      Accumulated or average volumetric outflow rate of solids
      real(c_real) BC_VOUT_s(DIMENSION_BC, DIM_M)
!
!                      Number of outflow rate values accumulated
      INTEGER          BC_OUT_N (DIMENSION_BC)
!
!                      Pressure drop specified for cyclic b.c. in X
      real(c_real) DELP_X
!
!                      Pressure drop specified for cyclic b.c. in Y
      real(c_real) DELP_Y
!
!                      Pressure drop specified for cyclic b.c. in Z
      real(c_real) DELP_Z
!
!                      Specified mass flux (e.g., g/cm^2.s) in the cyclic
!                      direction with specified pressure drop (only one
!                      direction is allowed).
      real(c_real) Flux_g
!
!                      Average gas velocity in X direction (for cyclic bc)
      real(c_real) U_g0
!
!                      Average gas velocity in Y direction (for cyclic bc)
      real(c_real) V_g0
!
!                      Average gas velocity in Z direction (for cyclic bc)
      real(c_real) W_g0
!
!                      Average solids velocity in X direction (for cyclic bc)
      real(c_real) U_s0 (DIM_M)
!
!                      Average solids velocity in Y direction (for cyclic bc)
      real(c_real) V_s0 (DIM_M)
!
!                      Average solids velocity in Z direction (for cyclic bc)
      real(c_real) W_s0 (DIM_M)
!
!                      IJK location where P_g is fixed for cyclic b.c's
      INTEGER          IJK_P_g(3)
!
!                      Coefficient in partial slip condition -- gas
      real(c_real) BC_hw_g (DIMENSION_BC)
!
!                      Coefficient in partial slip condition -- solids
      real(c_real) BC_hw_s (DIMENSION_BC, DIM_M)
!
!                      Wall velocity for partial slip condition -- gas
      real(c_real) BC_Uw_g (DIMENSION_BC)
!
!                      Wall velocity for partial slip condition -- gas
      real(c_real) BC_Vw_g (DIMENSION_BC)
!
!                      Wall velocity for partial slip condition -- gas
      real(c_real) BC_Ww_g (DIMENSION_BC)
!
!                      Wall velocity for partial slip condition -- solids
      real(c_real) BC_Uw_s (DIMENSION_BC, DIM_M)
!
!                      Wall velocity for partial slip condition -- solids
      real(c_real) BC_Vw_s (DIMENSION_BC, DIM_M)
!
!                      Wall velocity for partial slip condition -- solids
      real(c_real) BC_Ww_s (DIMENSION_BC, DIM_M)

!
!                      Coefficient in heat transfer boundary condition -- gas
      real(c_real) BC_hw_T_g (DIMENSION_BC)
!
!                      Coefficient in heat transfer boundary condition -- solids
      real(c_real) BC_hw_T_s (DIMENSION_BC, DIM_M)
!
!                      Wall temperature in heat transfer boundary  condition -- gas
      real(c_real) BC_Tw_g (DIMENSION_BC)
!
!                      Wall temperature in heat transfer boundary condition -- solids
      real(c_real) BC_Tw_s (DIMENSION_BC, DIM_M)
!
!                      Coefficient in heat transfer boundary condition -- gas
      real(c_real) BC_C_T_g (DIMENSION_BC)
!
!                      Coefficient in heat transfer boundary condition -- solids
      real(c_real) BC_C_T_s (DIMENSION_BC, DIM_M)

!
!                      Coefficient in mass transfer boundary condition -- gas
      real(c_real) BC_hw_X_g (DIMENSION_BC, DIM_N_g)
!
!                      Coefficient in mass transfer boundary condition -- solids
      real(c_real) BC_hw_X_s (DIMENSION_BC, DIM_M, DIM_N_s)
!
!                      Wall value in mass  transfer boundary  condition -- gas
      real(c_real) BC_Xw_g (DIMENSION_BC, DIM_N_g)
!
!                      Wall value in mass transfer boundary condition -- solids
      real(c_real) BC_Xw_s (DIMENSION_BC, DIM_M, DIM_N_s)
!
!                      Coefficient in mass transfer boundary condition -- gas
      real(c_real) BC_C_X_g (DIMENSION_BC, DIM_N_g)
!
!                      Coefficient in mass transfer boundary condition -- solids
      real(c_real) BC_C_X_s (DIMENSION_BC, DIM_M, DIM_N_s)
!

      LOGICAL:: CG_MI_CONVERTED_TO_PS(DIMENSION_BC)


! Flag to specify the constant number of particles per cell
! for the PIC solids
! Statistical weight of parcels will be calculated by the code
      INTEGER :: BC_PIC_MI_CONST_NPC(DIMENSION_BC, DIM_M)

! Flag to specify the constant statistical weight.
! for the PIC solids
! Number of computational particles/parcels will be calculated by the code
      real(c_real) :: BC_PIC_MI_CONST_STATWT(DIMENSION_BC, DIM_M)

      CONTAINS

         LOGICAL FUNCTION IS_CG(boundary_condition)
            implicit none
            INTEGER, intent(in) :: boundary_condition
            IS_CG = ((boundary_condition .eq. CG_PO) &
                 .or. (boundary_condition .eq. CG_MO) &
                 .or. (boundary_condition .eq. CG_NSW) &
                 .or. (boundary_condition .eq. CG_FSW) &
                 .or. (boundary_condition .eq. CG_PSW) &
                 .or. (boundary_condition .eq. CG_MI) &
                 )
         END FUNCTION IS_CG

         LOGICAL FUNCTION IS_NSW(boundary_condition)
            implicit none
            INTEGER, intent(in) :: boundary_condition
            IS_NSW = ((boundary_condition .eq. CG_NSW) &
                 )
         END FUNCTION IS_NSW

         LOGICAL FUNCTION IS_FSW(boundary_condition)
            implicit none
            INTEGER, intent(in) :: boundary_condition
            IS_FSW = ((boundary_condition .eq. CG_FSW) &
                 )
         END FUNCTION IS_FSW

         LOGICAL FUNCTION IS_PSW(boundary_condition)
            implicit none
            INTEGER, intent(in) :: boundary_condition
            IS_PSW = ((boundary_condition .eq. CG_PSW) )
         END FUNCTION IS_PSW

      END MODULE bc
