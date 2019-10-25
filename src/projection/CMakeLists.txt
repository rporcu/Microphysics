target_include_directories(mfixcore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    )

target_sources( mfixcore
    PRIVATE
 #   NodalProjection.H
 #   NodalProjection.cpp
    NodalProjector.H
    NodalProjector.cpp
    )
