add_sources ( mfix_eb.cpp get_eb_walls.f90 compute_normals.f90 )

# general EB geometry tools (projections onto edges, etc...)
add_sources ( eb_geometry.f90 )

# level-set sources
add_sources ( eb_levelset.cpp eb_levelset.H mfix_eb_F.H eb_levelset_F.f90 )
