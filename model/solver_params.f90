
      MODULE solver_params

      use bl_fort_module, only : c_real
      use iso_c_binding , only: c_int
 
      integer(c_int), parameter ::  pc_line = 0
      integer(c_int), parameter ::  pc_diag = 1
      integer(c_int), parameter ::  pc_none = 2

      integer(c_int), parameter ::  sweep_rsrs = 0
      integer(c_int), parameter ::  sweep_isis = 1
      integer(c_int), parameter ::  sweep_asas = 2

      END MODULE solver_params
