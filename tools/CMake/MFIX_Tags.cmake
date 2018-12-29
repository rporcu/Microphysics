#
# Macro to create tags
#
macro (add_tags_targets )

  if (( NOT DEFINED AMREX_INSTALL_DIR ) OR
      ( NOT PROJECT_SOURCE_DIR ) )
    message (AUTHOR_WARNING "Some paths are not defined")
    return ()
  endif ()


  set (SUPERBUILD_MODE ON)
  if (NOT EXISTS ${AMREX_INSTALL_DIR}/../sourcedir/Src )
    set (SUPERBUILD_MODE OFF)
  endif ()

  set (AMREX_SRC_DIR ${AMREX_INSTALL_DIR}/../sourcedir/Src/)

  if (ENABLE_PROJCC)
    set (MFIX_SRC_DIR  ${PROJECT_SOURCE_DIR}/src_cc/)
  else ()
    set (MFIX_SRC_DIR  ${PROJECT_SOURCE_DIR}/src_staggered/)
  endif ()

  set (MFIX_MODS_DIR ${PROJECT_SOURCE_DIR}/mods/)
  set (MFIX_DES_DIR  ${PROJECT_SOURCE_DIR}/src_des/)
  set (MFIX_EB_DIR   ${PROJECT_SOURCE_DIR}/src_eb/)

  set (TAGS_EXTRA_DIRS ${MFIX_MODS_DIR} ${MFIX_DES_DIR} ${MFIX_EB_DIR})

  # Add rule to generate TAGS
  # on macOS ctags-exuberant is just ctags
  find_program(CTAGS_EXU "ctags-exuberant")
  find_program(CTAGS_CMD "ctags")

  if(CTAGS_EXU)
    if (SUPERBUILD_MODE)
      add_custom_target ( tags
        COMMAND ctags-exuberant -R    --fortran-kinds=+i  ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMAND ctags-exuberant -R -e --fortran-kinds=+i  ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMENT "Generating tags files"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
      # Some file systems are case-insensitive so the above will fail:
      #  the first command would overwrite the second!
      add_custom_target ( ctags
        COMMAND ctags-exuberant -R    --fortran-kinds=+i ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMENT "Generating only ctags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
      add_custom_target( etags
        COMMAND ctags-exuberant -R -e --fortran-kinds=+i ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMENT "Generating only etags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    endif ()
    # Avoid "collisions" with AMReX => focus entirely on MFiX sources:
    add_custom_target ( mfix_ctags
      COMMAND ctags-exuberant -R    --fortran-kinds=+i ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
      COMMENT "Generating only ctags file for MFiX sources exclusively"
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    add_custom_target( mfix_etags
      COMMAND ctags-exuberant -R -e --fortran-kinds=+i ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
      COMMENT "Generating only etags file"
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
  elseif(CTAGS_CMD)
    if (SUPERBUILD_MODE)
      add_custom_target ( tags
        COMMAND ctags -R    --fortran-kinds=+i ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMAND ctags -R -e --fortran-kinds=+i ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMENT "Generating tags files"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
      # Some file systems are case-insensitive so the above will fail:
      #  the first command would overwrite the second!
      add_custom_target ( ctags
        COMMAND ctags -R    --fortran-kinds=+i ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMENT "Generating only ctags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
      add_custom_target( etags
        COMMAND ctags -R -e --fortran-kinds=+i ${AMREX_SRC_DIR} ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
        COMMENT "Generating only etags file"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    endif ()
    # Avoid "collisions" with AMReX => focus entirely on MFiX sources:
    add_custom_target ( mfix_ctags
      COMMAND ctags -R    --fortran-kinds=+i ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
      COMMENT "Generating only ctags file for MFiX sources exclusively"
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
    add_custom_target( mfix_etags
      COMMAND ctags-exuberant -R -e --fortran-kinds=+i ${MFIX_SRC_DIR} ${TAGS_EXTRA_DIRS}
      COMMENT "Generating only etags file"
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/ )
  endif()

  unset (SUPERBUILD_MODE)
  unset (MFIX_SRC_DIR)
  unset (AMREX_SRC_DIR)

endmacro ()
