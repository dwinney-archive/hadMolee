# Because the LoopTools library is build in fortran we need to compile the wrapper classes seperately
# and then link them to the core library later
# this way ROOT/cling doesnt get confused with mixed languages

# is there a better way to do this? maybe but i dont know it

# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "*.cpp")

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/LoopTools)
add_library( LOOPEXT SHARED ${SRC} )

set_property(TARGET LOOPEXT PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries( LOOPEXT ${LOOPLIB} )

install(TARGETS LOOPEXT
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )