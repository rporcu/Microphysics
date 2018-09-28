
  subroutine get_solver_params(eq_id, sweep_type, pc_type, max_it, tol) &
       bind(C, name="get_solver_params")

    use amrex_fort_module, only : rt => amrex_real
    use iso_c_binding , only: c_int, c_char

    use leqsol       , only: leq_sweep, leq_pc, leq_it, leq_tol
    use solver_params, only: sweep_rsrs, sweep_isis, sweep_asas
    use solver_params, only: pc_line, pc_diag, pc_none

    implicit none

    integer(c_int), intent(in   ) :: eq_id
    integer(c_int), intent(  out) :: sweep_type, pc_type, max_it
    real(rt)      , intent(  out) :: tol

    ! Default
    if (leq_sweep(eq_id) .eq. 'RSRS') then
       sweep_type = sweep_rsrs;

    elseif (leq_sweep(eq_id) .eq. 'ISIS') then
       sweep_type = sweep_isis;

    elseif (leq_sweep(eq_id) .eq. 'ASAS') then
       sweep_type = sweep_asas;

    else
       print *,'Dont know this leq_sweep flag'
       stop
    end if

    ! Default
    if (leq_pc(eq_id) .eq. 'LINE') then
       pc_type = pc_line;

    elseif (leq_pc(eq_id) .eq. 'DIAG') then
       pc_type = pc_diag;

    elseif (leq_pc(eq_id) .eq. 'NONE') then
       pc_type = pc_none;

    else
       print *,'Dont know this leq_sweep flag'
       stop
    end if

    max_it = leq_it(eq_id);
    tol    = leq_tol(eq_id);

  end subroutine get_solver_params
