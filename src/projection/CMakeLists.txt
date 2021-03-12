target_include_directories(mfixcore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    )

target_sources( mfixcore
    PRIVATE
    mfix_apply_nodal_projection.cpp
    mfix_proj_redistribution.cpp
    mfix_constraints.cpp
    )
