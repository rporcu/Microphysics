//
//  Redistribution.H
//
#ifndef redistribution_H
#define redistribution_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>

namespace redistribution {

    void redistribute_eb (amrex::Box const& bx, int ncomp, int scomp,
                          amrex::Array4<amrex::Real>       const& dUdt_out,
                          amrex::Array4<amrex::Real const> const& dUdt_in,
                          amrex::Array4<amrex::Real const> const& U_in,
                          amrex::Array4<amrex::Real> const& scratch,
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                       amrex::Array4<amrex::Real const> const& vmac,
                                       amrex::Array4<amrex::Real const> const& wmac),
                          amrex::Array4<amrex::EBCellFlag const> const& flag,
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                       amrex::Array4<amrex::Real const> const& apy,
                                       amrex::Array4<amrex::Real const> const& apz),
                          amrex::Array4<amrex::Real const> const& vfrac,
                          AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                       amrex::Array4<amrex::Real const> const& fcy,
                                       amrex::Array4<amrex::Real const> const& fcz),
                          amrex::Array4<amrex::Real const> const& ccent,
                          amrex::Geometry& geom, amrex::Real dt, std::string advection_type);

    void flux_redistribute_eb ( amrex::Box const& bx, int ncomp, int scomp,
                                amrex::Array4<amrex::Real>       const& dUdt_out,
                                amrex::Array4<amrex::Real const> const& dUdt_in,
                                amrex::Array4<amrex::Real> const& scratch,
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
                                amrex::Array4<amrex::Real const> const& vfrac,
                                amrex::Geometry& geom);
    void state_redistribute_eb ( amrex::Box const& bx, int ncomp, int scomp,
                                 amrex::Array4<amrex::Real>       const& dUdt_out,
                                 amrex::Array4<amrex::Real const> const& dUdt_in,
                                 amrex::Array4<amrex::EBCellFlag const> const& flag,
                                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                              amrex::Array4<amrex::Real const> const& apy,
                                              amrex::Array4<amrex::Real const> const& apz),
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                              amrex::Array4<amrex::Real const> const& fcy,
                                              amrex::Array4<amrex::Real const> const& fcz),
                                 amrex::Array4<amrex::Real const> const& ccent,
                                 amrex::Geometry& geom);

    void make_itracker             ( amrex::Box const& bx,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                  amrex::Array4<amrex::Real const> const& apy,
                                                  amrex::Array4<amrex::Real const> const& apz),
                                     amrex::Array4<amrex::Real const> const& vfrac,
                                     amrex::Array4<int> const& itracker,
                                     amrex::Geometry& geom);

    void merge_redistribute_update ( amrex::Box const& bx, int ncomp, int scomp,
                                     amrex::Array4<amrex::Real>       const& dUdt_out,
                                     amrex::Array4<amrex::Real const> const& dUdt_in,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                  amrex::Array4<amrex::Real const> const& apy,
                                                  amrex::Array4<amrex::Real const> const& apz),
                                     amrex::Array4<amrex::Real const> const& vfrac,
                                     amrex::Array4<int> const& itracker,
                                     amrex::Geometry& geom);

    void merge_redistribute_full   ( amrex::Box const& bx, int ncomp, int scomp,
                                     amrex::Array4<amrex::Real>       const& dUdt_out,
                                     amrex::Array4<amrex::Real const> const& dUdt_in,
                                     amrex::Array4<amrex::Real const> const& U_in,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                                  amrex::Array4<amrex::Real const> const& vmac,
                                                  amrex::Array4<amrex::Real const> const& wmac),
                                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                                  amrex::Array4<amrex::Real const> const& apy,
                                                  amrex::Array4<amrex::Real const> const& apz),
                                     amrex::Array4<amrex::Real const> const& vfrac,
                                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                                  amrex::Array4<amrex::Real const> const& fcy,
                                                  amrex::Array4<amrex::Real const> const& fcz),
                                     amrex::Array4<amrex::Real const> const& ccent,
                                     amrex::Geometry& geom, amrex::Real dt);











  void apply_eb_redistribution ( const amrex::Box& bx,
                                 amrex::Array4<amrex::Real>       const& div,
                                 amrex::Array4<amrex::Real const> const& divc,
                                 amrex::Array4<amrex::Real const> const& weights,
                                 amrex:: MFIter const& mfi,
                                 const int icomp,
                                 const int ncomp,
                                 amrex::Array4<amrex::EBCellFlag const> const& flag,
                                 amrex::Array4<amrex::Real const> const& vfrac,
                                 const amrex::Geometry & geom);









} // namespace redistribution

#endif /* Redistribution_H */
