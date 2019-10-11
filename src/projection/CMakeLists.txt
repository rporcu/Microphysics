target_include_directories(mfixcore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    )

target_sources( mfixcore
    PRIVATE
    MFIX_NodalProjection.H
    MFIX_NodalProjection.cpp
    )
