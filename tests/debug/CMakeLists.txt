set ( DEBUG_TESTS FLD05-xy FLD05-xz FLD05-yz )

set (MFIX_CMD ${CMAKE_BINARY_DIR}/mfix )

foreach ( THIS_TEST ${DEBUG_TESTS} )

   # string ( REPLACE "/" "_" DEBUG_TEST_NAME ${THIS_TEST} )
   set ( DEBUG_TEST_NAME "run_${THIS_TEST}" )
   set ( DEST_DIR ${CMAKE_BINARY_DIR}/debug/${THIS_TEST} )
   set ( SOURCE_DIR ${PROJECT_SOURCE_DIR}/debug/${THIS_TEST} )
   
   # Copy file to binary dir
   file ( COPY ${SOURCE_DIR} DESTINATION ${CMAKE_BINARY_DIR}/debug/ )

   add_custom_target ( ${DEBUG_TEST_NAME}
      DEPENDS mfix
      COMMAND [ -d ${DEST_DIR} ] || mkdir ${DEST_DIR}
      COMMAND cp -R --update ${SOURCE_DIR}/* .
      COMMAND ${MFIX_CMD} inputs
      WORKING_DIRECTORY  ${DEST_DIR} )
   
endforeach ()

