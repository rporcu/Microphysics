################################################################
#                                                              #
#      Automake input file                                     #
#                                                              #
################################################################

# Automake defaults to F77LINK.  Using FC to link will work whether FC is an MPI or non-MPI compiler.
F77LD = $(FC)

PPFCCOMPILE = $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(FPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(FCFLAGS_f)

.F.o:
	$(F77COMPILE) -c -o $@ $<

.f.o:
	$(PPFCCOMPILE) -c -o $@ $<

FC_MODEXT=mod

# directory with all the *.mod files
MODDIR=$(top_builddir)/model
MODOUTDIR=$(MODDIR)
CLEANFILES=$(MODOUTDIR)/*.$(FC_MODEXT) $(libmfix_src:%.f=%.d)
DISTCLEANFILES = build-aux/Makefile.usr

AM_FCFLAGS = $(FC_MODINC)$(top_srcdir)/model/include
AM_FCFLAGS += $(FC_MODINC)$(MODDIR)
AM_FCFLAGS += $(FC_MODOUT)$(MODOUTDIR)
AM_LDFLAGS =


SUFFIXES = .d

.DEFAULT_GOAL := mfix$(EXEEXT)

################################################################
#                                                              #
#           mfix executable                                    #
#                                                              #
################################################################
bin_PROGRAMS = mfix$(EXEEXT)

mfix_LDADD = libmfix.a
noinst_LIBRARIES = libmfix.a

# mfix.f needs to be in mfix_SOURCES (not just libmfix_a_SOURCES) for gfortran 4.4
mfix_SOURCES = model/mfix.f
libmfix_a_SOURCES =


################################################################
#                                                              #
#      EXTRA_DIST:  Non-source and non-binary files that       #
#                   shoud be included in "make dist" tarball   #
#                                                              #
################################################################

EXTRA_DIST =  Bug_Fixes CHANGES DISCLAIMER NOTICE benchmarks build-aux config configure_mfix doc
EXTRA_DIST += legacy_tests mfix_user_guide.pdf readme regression_tests tests tools tutorials visualization_tools
EXTRA_DIST += model/make_mfix
EXTRA_DIST += model/make_mfix.old

################################################################
#                                                              #
#  Use script "fortran-depcomp" to find dependencies among     #
#  source files.                                               #
#                                                              #
################################################################

.f.d:
	@$(MKDIR_P) $(@D)
	SRCDIR="$(srcdir)" OBJEXT="$(OBJEXT)" MODEXT="$(FC_MODEXT)" MODDIR="$(MODDIR)" \
	SED="$(SED)" GREP="$(GREP)" $(SHELL) $(top_srcdir)/build-aux/fortran-depcomp $^ > $@

# Use generated module dependency rules

@make_ifnotclean@
$(eval -include $(libmfix_src:%.f=%.d) )
@make_endif@

################################################################
#      TAU support
################################################################
TAU_F90=@TAU_F90@
TAU_MAKEFILE=@TAU_MAKEFILE@
export TAU_F90
export TAU_MAKEFILE


################################################################
#                                                              #
#      MKL support (only with Intel Fortran compiler)          #
#                                                              #
################################################################

if MKL
AM_FCFLAGS += -mkl
AM_LDFLAGS += -mkl
else
libmfix_a_SOURCES += model/solvers/BLAS.f
libmfix_a_SOURCES += model/solvers/DGTSV.f
endif

################################################################
#                                                              #
#      OpenMP and MPI support                                  #
#                                                              #
#      --enable-smp and --enable-dmp options in configure.ac   #
#                                                              #
################################################################

if OPENMP
AM_FCFLAGS += $(OPENMP_FCFLAGS)
AM_LDFLAGS += $(OPENMP_FCFLAGS)
endif

if MPI
AM_FCFLAGS += $(FC_DEFINE)MPI
AM_LDFLAGS += $(FC_DEFINE)MPI
endif

################################################################
#      libmfix_a_SOURCES                                       #
#                                                              #
#      The majority of source files for mfix                   #
#      are listed here.                                        #
################################################################

libmfix_a_SOURCES += \
	model/include/boundary_conditions.inc \
	model/include/cartesian_grid_namelist.inc \
	model/include/chem_equations.inc \
	model/include/desnamelist.inc \
	model/include/dmp_batch_control.inc \
	model/include/fun_avg.inc \
	model/include/functions.inc \
	model/include/gas_phase.inc \
	model/include/geometry.inc \
	model/include/initial_conditions.inc \
	model/include/legacy.inc \
	model/include/numerical_params.inc \
	model/include/output_control.inc \
	model/include/physical_params.inc \
	model/include/point_sources.inc \
	model/include/qmomknamelist.inc \
	model/include/run_control.inc \
	model/include/solids_phase.inc \
	model/include/species.inc \
	model/include/species_indices.inc \
	model/include/tfm_solids.inc \
	model/include/usr_hooks.inc \
	model/include/usrnlst.inc \
	model/include/xforms.inc

libmfix_src = \
	model/accum_resid.f \
	model/adjust_a_g.f \
	model/adjust_dt_mod.f \
	model/adjust_leq.f \
	model/allocate_arrays.f \
	model/bc_mod.f \
	model/bodyforce_mod.f \
	model/boundfunijk_mod.f \
	model/calc_cell.f \
	model/calc_coeff.f \
	model/calc_d_e.f \
	model/calc_d_n.f \
	model/calc_d_t.f \
	model/calc_mflux.f \
	model/calc_mu_g.f \
	model/calc_outflow.f \
	model/calc_resid.f \
	model/calc_trd_g.f \
	model/cartesian_grid/CG_set_bc0.f \
	model/cartesian_grid/CG_set_outflow.f \
	model/cartesian_grid/CG_source_u_g.f \
	model/cartesian_grid/CG_source_v_g.f \
	model/cartesian_grid/CG_source_w_g.f \
	model/cartesian_grid/allocate_cut_cell_arrays.f \
	model/cartesian_grid/allocate_dummy_cut_cell_arrays.f \
	model/cartesian_grid/calc_vort_out.f \
	model/cartesian_grid/cartesian_grid_init_namelist.f \
	model/cartesian_grid/check_data_cartesian.f \
	model/cartesian_grid/cut_cell_preprocessing.f \
	model/cartesian_grid/cutcell_mod.f \
	model/cartesian_grid/dashboard_mod.f \
	model/cartesian_grid/deallocate_cut_cell_arrays.f \
	model/cartesian_grid/define_quadrics.f \
	model/cartesian_grid/dmp_cartesian.f \
	model/cartesian_grid/eval_usr_fct.f \
	model/cartesian_grid/get_alpha.f \
	model/cartesian_grid/get_connectivity.f \
	model/cartesian_grid/get_cut_cell_flags.f \
	model/cartesian_grid/get_cut_cell_volume_area.f \
	model/cartesian_grid/get_delh.f \
	model/cartesian_grid/get_master.f \
	model/cartesian_grid/get_poly_data.f \
	model/cartesian_grid/get_stl_data.f \
	model/cartesian_grid/polygon_mod.f \
	model/cartesian_grid/progress_bar_mod.f \
	model/cartesian_grid/quadric_mod.f \
	model/cartesian_grid/set_Odxyz.f \
	model/cartesian_grid/stl_mod.f \
	model/cartesian_grid/update_dashboard.f \
	model/cartesian_grid/vtk_mod.f \
	model/cartesian_grid/vtk_out.f \
	model/cartesian_grid/write_progress_bar.f \
	model/check_ab_m.f \
	model/check_bqend.f \
	model/check_convergence.f \
	model/check_data/check_axis.f \
	model/check_data/check_bc_dem.f \
	model/check_data/check_bc_geometry.f \
	model/check_data/check_bc_inflow.f \
	model/check_data/check_bc_outflow.f \
	model/check_data/check_bc_walls.f \
	model/check_data/check_boundary_conditions.f \
	model/check_data/check_dmp_prereqs.f \
	model/check_data/check_gas_phase.f \
	model/check_data/check_geometry.f \
	model/check_data/check_geometry_prereqs.f \
	model/check_data/check_ic_common_discrete.f \
	model/check_data/check_ic_dem.f \
	model/check_data/check_initial_conditions.f \
	model/check_data/check_numerics.f \
	model/check_data/check_output_control.f \
	model/check_data/check_point_sources.f \
	model/check_data/check_run_control.f \
	model/check_data/check_solids_common_all.f \
	model/check_data/check_solids_common_discrete.f \
	model/check_data/check_solids_dem.f \
	model/check_data/check_solids_model_prereqs.f \
	model/check_data/check_solids_phases.f \
	model/check_data_20.f \
	model/check_data_30.f \
	model/check_plane.f \
	model/constant_mod.f \
	model/conv_dif_u_g.f \
	model/conv_dif_v_g.f \
	model/conv_dif_w_g.f \
	model/conv_pp_g.f \
	model/conv_rop.f \
	model/corner.f \
	model/correct_0.f \
	model/dbg_mod.f \
	model/deprecated.f \
	model/des/calc_collision_wall_mod.f \
	model/des/calc_drag_des.f \
	model/des/calc_epg_des.f \
	model/des/calc_force_dem.f \
	model/des/calc_grad_des.f \
	model/des/calc_interp_weights.f \
	model/des/calc_pg_grad.f \
	model/des/cfassign.f \
	model/des/cffctowall.f \
	model/des/cfnewvalues.f \
	model/des/cfrelvel.f \
	model/des/cfslide.f \
	model/des/cfupdateold.f \
	model/des/check_cell_movement.f \
	model/des/comp_mean_fields.f \
	model/des/comp_mean_fields0.f \
	model/des/comp_mean_fields1.f \
	model/des/des_allocate_mod.f \
	model/des/des_bc_mod.f \
	model/des/des_granular_temperature.f \
	model/des/des_init_arrays.f \
	model/des/des_init_namelist.f \
	model/des/des_time_march.f \
	model/des/desgrid_mod.f \
	model/des/desmpi_mod.f \
	model/des/desmpi_wrapper_mod.f \
	model/des/discretelement_mod.f \
	model/des/drag_gp_des.f \
	model/des/drag_gs_des0.f \
	model/des/drag_gs_des1.f \
	model/des/gas_drag.f \
	model/des/generate_particles_mod.f \
	model/des/interpolation_mod.f \
	model/des/layout_mi_dem.f \
	model/des/make_arrays_des.f \
	model/des/mass_inflow_dem.f \
	model/des/mass_outflow_dem.f \
	model/des/mpi_comm_des_mod.f \
	model/des/mpi_funs_des_mod.f \
	model/des/mpi_init_des_mod.f \
	model/des/mpi_node_des_mod.f \
	model/des/mpi_pack_des_mod.f \
	model/des/mpi_unpack_des_mod.f \
	model/des/neighbour.f \
	model/des/nsquare.f \
	model/des/particle_filter_mod.f \
	model/des/particles_in_cell.f \
	model/des/randomno_mod.f \
	model/des/read_particle_input.f \
	model/des/read_res0_des.f \
	model/des/read_res1_des_mod.f \
	model/des/sendrecvnode_mod.f \
	model/des/set_bc_dem.f \
	model/des/set_bc_dem_mi.f \
	model/des/set_bc_dem_mo.f \
	model/des/set_filter_des.f \
	model/des/set_geometry_des.f \
	model/des/set_phase_index.f \
	model/des/solid_drag.f \
	model/des/stl_dbg_des_mod.f \
	model/des/stl_functions_des_mod.f \
	model/des/stl_preproc_des_mod.f \
	model/usr/usr0_des.f \
	model/usr/usr1_des.f \
	model/usr/usr2_des.f \
	model/usr/usr3_des.f \
	model/des/vtp_mod.f \
	model/des/write_des_data.f \
	model/des/write_res0_des.f \
	model/des/write_res1_des_mod.f \
	model/solvers/dgtsl.f \
	model/discretization_mod.f \
	model/display_resid.f \
	model/dmp_modules/compar_mod.f \
	model/dmp_modules/dbg_util_mod.f \
	model/dmp_modules/debug_mod.f \
	model/dmp_modules/gridmap_mod.f \
	model/dmp_modules/mpi_mod.f \
	model/dmp_modules/mpi_utility_mod.f \
	model/dmp_modules/parallel_mpi_mod.f \
	model/dmp_modules/sendrecv3_mod.f \
	model/dmp_modules/sendrecv_mod.f \
	model/drag_mod.f \
	model/eos_mod.f \
	model/error_manager_mod.f \
	model/exit.f \
	model/fldvar_mod.f \
	model/flow_to_vel.f \
	model/fun_avg_mod.f \
	model/functions_mod.f \
	model/io/funits_mod.f \
	model/geometry_mod.f \
	model/get_bc_area.f \
	model/get_data.f \
	model/get_ps.f \
	model/ic_mod.f \
	model/io/in_binary_512_mod.f \
	model/io/in_binary_512i_mod.f \
	model/indices_mod.f \
	model/init_fvars.f \
	model/init_namelist.f \
	model/iterate.f \
	model/solvers/leq_bicgs.f \
	model/solvers/leq_bicgst.f \
	model/solvers/leq_cg.f \
	model/solvers/leq_gmres.f \
	model/solvers/leq_sor.f \
	model/solvers/leqsol_mod.f \
	model/location_check.f \
	model/machine_mod.f \
	model/matrix_mod.f \
	model/mfix.f \
	model/mflux_mod.f \
	model/mod_bc_i.f \
	model/mod_bc_j.f \
	model/mod_bc_k.f \
	model/io/open_file.f \
	model/io/open_files.f \
	model/io/out_array.f \
	model/io/out_array_c.f \
	model/io/out_array_k.f \
	model/io/out_array_kc.f \
	model/io/out_bin_512.f \
	model/io/out_bin_512i.f \
	model/io/out_bin_512r.f \
	model/io/out_bin_r.f \
	model/io/output_manager.f \
	model/output_mod.f \
	model/param1_mod.f \
	model/param_mod.f \
	model/parse_line.f \
	model/parse_resid_string.f \
	model/physical_prop.f \
	model/physprop_mod.f \
	model/ps_mod.f \
	model/read_namelist.f \
	model/io/read_res1.f \
	model/remove_comment.f \
	model/reset_new.f \
	model/residual_mod.f \
	model/run_mod.f \
	model/scales_mod.f \
	model/set_bc0.f \
	model/set_bc1.f \
	model/set_bc_flow.f \
	model/set_constants.f \
	model/set_constprop.f \
	model/set_flags.f \
	model/set_fluidbed_p.f \
	model/set_geometry.f \
	model/set_geometry1.f \
	model/set_ic.f \
	model/set_icbc_flags.f \
	model/set_increments.f \
	model/set_index1.f \
	model/set_index1a.f \
	model/set_max2.f \
	model/set_outflow.f \
	model/set_param.f \
	model/set_ps.f \
	model/set_ro_g.f \
	model/set_wall_bc.f \
	model/shift_dxyz.f \
	model/solvers/solve_lin_eq.f \
	model/solve_pp_g.f \
	model/solve_vel_star.f \
	model/source_pp_g.f \
	model/source_u_g.f \
	model/source_v_g.f \
	model/source_w_g.f \
	model/tau_u_g.f \
	model/tau_v_g.f \
	model/tau_w_g.f \
	model/solvers/test_lin_eq.f \
	model/time_cpu_mod.f \
	model/time_march.f \
	model/tmp_array1_mod.f \
	model/tmp_array_mod.f \
	model/toleranc_mod.f \
	model/under_relax.f \
	model/update_old.f \
	model/ur_facs_mod.f \
	model/usr/usr0.f \
	model/usr/usr1.f \
	model/usr/usr2.f \
	model/usr/usr3.f \
	model/usr/usr_drag.f \
	model/usr/usr_init_namelist.f \
	model/usr/usr_mod.f \
	model/usr/usr_physical_prop.f \
	model/usr/usr_write_out0.f \
	model/usr/usr_write_out1.f \
	model/utilities_mod.f \
	model/vavg_u_g.f \
	model/vavg_v_g.f \
	model/vavg_w_g.f \
	model/vf_gs.f \
	model/io/write_ab_m.f \
	model/io/write_ab_m_var.f \
	model/io/write_error.f \
	model/io/write_header.f \
	model/io/write_out0.f \
	model/io/write_out1.f \
	model/io/write_out3.f \
	model/io/write_res0.f \
	model/io/write_res1.f \
	model/io/write_spx0.f \
	model/io/write_spx1.f \
	model/io/write_table.f \
	model/usr/write_usr0.f \
	model/usr/write_usr1.f \
	model/solvers/xerbla.f \
	model/xsi_array_mod.f \
	model/xsi_mod.f \
	model/zero_norm_vel.f

libmfix_a_SOURCES += $(libmfix_src)


cases := $(shell find . -name mfix.dat)

allcases:
	for case in $(cases); do \
		cp build-aux/Makefile.usr $$(dirname $$case)/Makefile || exit; \
		$(MAKE) -C $$(dirname $$case) mfix$(EXEEXT) || exit; \
	done
