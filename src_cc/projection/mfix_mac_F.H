#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include "AMReX_RealVect.H"

#ifdef __cplusplus
extern "C"
{
#endif


    void compute_ugradu_mac (
	const int* lo, const int* hi,
	amrex::Real* rhs, const int* rlo, const int* rhi, 
	const amrex::Real* vel   ,  const int* vello, const int* velhi,
	const amrex::Real* ep_g   ,  const int* elo, const int* ehi,
	amrex::Real* u, const int* ulo, const int* uhi,
	amrex::Real* v, const int* vlo, const int* vhi,	
	amrex::Real* w, const int* wlo, const int* whi,
	const amrex::Real* xslopes, 
	const amrex::Real* yslopes, 
	const amrex::Real* zslopes, const int* slo, const int* shi, 
	const int* domlo, const int* domhi,
	const int* bc_ilo_type, const int* bc_ihi_type, 
	const int* bc_jlo_type, const int* bc_jhi_type,
	const int* bc_klo_type, const int* bc_khi_type,
	const amrex::Real* dx, const int* ng, const int* ugradu_type
        );

    void compute_mac_diveu (
	const int* lo, const int* hi,
	amrex::Real* diveu,  const int* slo, const int* shi,
	const amrex::Real* ep_g,
	const amrex::Real* u_g, const int* ulo, const int* uhi,
	const amrex::Real* v_g, const int* vlo, const int* vhi,
	const amrex::Real* w_g, const int* wlo, const int* whi,
	const amrex::Real* dx
        );

    void project_mac_velocity (
	const int* lo, const int* hi,
	amrex::Real* u_i, const int* ulo, const int* uhi,
	const amrex::Real* dphidxi, const int* glo, const int* ghi,
	const amrex::Real* ro, const int* slo, const int* shi,
	const amrex::Real* c, const int* dir
        );

    void set_mac_velocity_bcs( 
        const int* slo, const int* shi,
        amrex::Real* u_g,  const int* ulo, const int* uhi,
        amrex::Real* v_g,  const int* vlo, const int* vhi,
        amrex::Real* w_g,  const int* wlo, const int* whi,
        const int* bc_ilo, const int* bc_ihi,
        const int* bc_jlo, const int* bc_jhi,
        const int* bc_klo, const int* bc_khi,
        const int* domlo , const int* domhi,
        const int* ng
        );

    void compute_bcoeff_mac (
	const int* lo, const int* hi,
	amrex::Real* bcoeff, const int* alo, const int* ahi,
	const amrex::Real* u,     const int* ulo, const int* uhi,
	const amrex::Real* ro_g, const int* slo, const int* shi,
	const amrex::Real* ep_g,
	const int* dir
        );

#ifdef __cplusplus
			       }
#endif
