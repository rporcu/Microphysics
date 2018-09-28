module avg_functions

   use amrex_fort_module, only : rt => amrex_real

   contains

  ! Arithmetic average
   real(rt) function avg(Xm, Xp)
      implicit none
      double precision, intent(in) :: Xp, Xm
      avg = 0.5d0 *(Xm + Xp)
   end function avg

  ! Harmonic average
   real(rt) function avg_h(Xm, Xp)
      use param, only: small_number
      implicit none
      double precision, intent(in) :: Xp, Xm
      avg_h = Xm * Xp / max(small_number, 0.5d0*(Xm + Xp) )
   end function avg_h

end module avg_functions
