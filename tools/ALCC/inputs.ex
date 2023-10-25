# Parameters for app

# MFIX-Exa plot file to process
#plotfile = /nfs/home/2/jmusser/projects/ALCC/plt18626
#plotfile = /nfs/home/2/jmusser/projects/ALCC/ORG-PLOTS/plt43668
#plotfile = /nfs/home/2/jmusser/projects/ALCC/plt70192
plotfile = /nfs/home/2/jmusser/projects/ALCC/phi0.001_ar1.8/plt00000

# Sets the max index for looping through multiple plot files. The
# name (and path) are defined by the plotfile keyword. If you
# want to loop over multiple files, set the following to the largest
# plot file index.
#
#plotfile.max_index = 500 # int [-1]

#==============================================================================
# Defines the mean particle diameter which may be needed for
# various calculations.
#
#  * If the following is left undefined, then it is computed
#    using particle radii values in the plot file.
#  * A provided value overrides a computed mean diameter.
#
# particle.mean_diameter = Real [computed / undefined]

#==============================================================================
# Convert native plot file to HDF5. This requires that post-mfix was built
# with HDF5 support.
#
# This requires the app
#hdf5 = true # bool [false]
#
# Option to use data compression.
# See notes in  /amrex/Tests/HDF5Benchmark/GNUMakefile for details
#hdf5.compression = None@0 # string [None@0], ZFP_ACCURACY@0.001
#
# By default, HDF5 files are written to the same directory that contains
# the plot file. By specifying a path, the output can be redirected to
# a different directory.
#hdf5.path = /nfs/home/2/jmusser/tmp # string [./]

#==============================================================================
# Compute global avareges of known quanties
#
# It is assumed that the fluid velocity components,
# pressure and pressure gradient are saved in the
# plot file in addition to particle radi and
# velocity components.
#
# terms (case sensitive)
#..............................................................................
#
# alpha_f    -- fluid volume fraction
# Uf, Vf, Wf -- fluid velocity components
# Pf         -- fluid pressure
# sigma_##   -- component of stress tensor (11 -> 33)
#
# alpha_p    -- solids volume fraction
# Up, Vp, Wp -- particle velocity components (Eulerian)
# sigPp_##   -- granular viscosity terms (11 -> 33)
#
# operations (case sensitive)
#..............................................................................
#
# grad_x()   -- d/dx of term in parenthesis
# grad_y()   -- d/dy of term in parenthesis
# grad_z()   -- d/dz of term in parenthesis
#
# A*B        -- Product of terms A and B
#
#
averages = none
#averages = alpha_p Up Vp Wp \
#           alpha_f Uf Vf Wf Pf \
#                               \
#           Pf*grad_x(Uf) grad_x(alpha_f*Uf) \
#           Pf*grad_y(Vf) grad_y(alpha_f*Vf) \
#           Pf*grad_z(Wf) grad_z(alpha_f*Wf) \
#                                            \
#           sigma_11*grad_x(Uf) sigma_11 grad_x(Uf) \
#           sigma_12*grad_y(Uf) sigma_12 grad_y(Uf) \
#           sigma_13*grad_z(Uf) sigma_13 grad_z(Uf) \
#                                                   \
#           sigma_21*grad_x(Vf) sigma_21 grad_x(Vf) \
#           sigma_22*grad_y(Vf) sigma_22 grad_y(Vf) \
#           sigma_23*grad_z(Vf) sigma_23 grad_z(Vf) \
#                                                   \
#           sigma_31*grad_x(Wf) sigma_31 grad_x(Wf) \
#           sigma_32*grad_y(Wf) sigma_32 grad_y(Wf) \
#           sigma_33*grad_z(Wf) sigma_33 grad_z(Wf) \
#                                                   \
#           alpha_p*Up*Uf  alpha_p*Up  alpha_p*Uf alpha_p*Uf*Uf  alpha_f*Uf \
#           alpha_p*Vp*Vf  alpha_p*Vp  alpha_p*Vf alpha_p*Vf*Vf  alpha_f*Vf \
#           alpha_p*Wp*Wf  alpha_p*Wp  alpha_p*Wf alpha_p*Wf*Wf  alpha_f*Wf \
#                                                                           \
#           grad_x(Pf) Uf*grad_x(Pf) \
#           grad_y(Pf) Vf*grad_y(Pf) \
#           grad_z(Pf) Wf*grad_z(Pf) \
#                                    \
#           alpha_p*Uf*grad_x(sigma_11)  grad_x(sigma_11) \
#           alpha_p*Uf*grad_x(sigma_12)  grad_x(sigma_12) \
#           alpha_p*Uf*grad_x(sigma_13)  grad_x(sigma_13) \
#                                                         \
#           alpha_p*Vf*grad_y(sigma_21)  grad_y(sigma_21) \
#           alpha_p*Vf*grad_y(sigma_22)  grad_y(sigma_22) \
#           alpha_p*Vf*grad_y(sigma_23)  grad_y(sigma_23) \
#                                                         \
#           alpha_p*2f*grad_z(sigma_31)  grad_z(sigma_31) \
#           alpha_p*2f*grad_z(sigma_32)  grad_z(sigma_32) \
#	   alpha_p*2f*grad_z(sigma_33)  grad_z(sigma_33) \
#                                                         \
#           alpha_p*Theta alpha_p*Theta*grad_x(Uf) \
#                         alpha_p*Theta*grad_y(Vf) \
#                         alpha_p*Theta*grad_z(Wf) \
#                                                  \
#           alpha_p*sigPp_11*grad_x(Uf)  sigPp_11 \
#           alpha_p*sigPp_12*grad_y(Uf)  sigPp_12 \
#           alpha_p*sigPp_13*grad_z(Uf)  sigPp_13 \
#                                                 \
#           alpha_p*sigPp_21*grad_x(Vf)  sigPp_21 \
#           alpha_p*sigPp_22*grad_y(Vf)  sigPp_22 \
#           alpha_p*sigPp_23*grad_z(Vf)  sigPp_23 \
#                                                 \
#           alpha_p*sigPp_31*grad_x(Wf)  sigPp_31 \
#           alpha_p*sigPp_32*grad_y(Wf)  sigPp_32 \
#           alpha_p*sigPp_33*grad_z(Wf)  sigPp_33 \
#                                                 \
#           Up*grad_x(Pf) alpha_p*grad_x(Pf) \
#           Vp*grad_y(Pf) alpha_p*grad_y(Pf) \
#           Wp*grad_z(Pf) alpha_p*grad_z(Pf) \
#                                            \
#           alpha_p*Up*grad_x(sigma_11) \
#           alpha_p*Vp*grad_y(sigma_22) \
#           alpha_p*Wp*grad_z(sigma_33)


# Options for filtering deposited particle data:
# ------------------------------------------------------------//
# The filer can be defined as constant (no smoothing),
# constant or variable.

filter.type = none # string [none], constant, variable

# Constant filter size .........................................
#
# The filter width can be explicitly defined. The provided value
# specifies the diffusion coefficient to be used as-is. This
# input overrides a defined width (filter.width)
#
# filter.diff_coeff = Real [undefined]
#
# The filter size can also be defined as a function of the mean
# particle diameter. (See particle.mean_diameter)
#
# filter.width = 7.0 # Real [7.0]
#
# Variable filter size .........................................
#
# A local filter size is computed from the solids volume fraction
#
# filter.sample_size = 10.0 # Real [10.0]


#==============================================================================
# Settings for diffusion solver
# ------------------------------------------------------------//
diffusion.verbose        = 0 # int [0]
diffusion.bottom_verbose = 0 # int [0]

diffusion.max_iter        = 500 # int [100]
diffusion.bottom_max_iter = 100 # int [100]

diffusion.max_fmg_iter   = 0 # int [0]
diffusion.linop_maxorder = 2 # int [2]

diffusion.agglomeration  = true  # bool [true]
diffusion.agg_grid_size  = -1

diffusion.consolidation  = true  # bool [true]
diffusion.semicoarsening = false # bool [false]

diffusion.max_coarsening_level     = 32 # int [32]
diffusion.max_semicoarsening_level = 32 # int [32]

diffusion.rtol = 1.0e-4  # Real [1.0e-11]
diffusion.atol = 1.0e-9  # Real [1.0e-14]

diffusion.bottom_solver = "bicgstab" # string ["bicgstab"]
