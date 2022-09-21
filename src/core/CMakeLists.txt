# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "*.cpp")

if (BUILD_CLING)
    ROOT_GENERATE_DICTIONARY(G__EE3BODY
            reaction_kinematics.hpp
            amplitude.hpp
            amplitude_sum.hpp
            breit_wigner.hpp
            triangle.hpp
            hadronic_molecule.hpp
            LINKDEF EE3BODYLinkDef.h
            )
    add_library( EE3BODY SHARED ${SRC} G__EE3BODY.cxx )
else()
    add_library( EE3BODY SHARED ${SRC} )
endif()

target_link_libraries( EE3BODY ROOT::Core ROOT::MathMore )
target_link_libraries( EE3BODY ${CUBATURELIB} )

set_property(TARGET EE3BODY PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries( EE3BODY ${LOOPLIB} )

install(TARGETS EE3BODY
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

if (BUILD_CLING)
    install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libEE3BODY_rdict.pcm    
            DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()