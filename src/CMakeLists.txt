# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "core/*.cpp")

include_directories("core")

if (BUILD_DICTIONARY)
    ROOT_GENERATE_DICTIONARY(G__HADMOLEE
            kinematics.hpp
            amplitude.hpp
            lineshape.hpp
            breit_wigner.hpp
            triangle.hpp
            box.hpp
            molecule.hpp
            fitter.hpp
            LINKDEF HADMOLEELinkDef.h
            )
    add_library( HADMOLEE SHARED ${SRC} G__HADMOLEE.cxx )
else()
    add_library( HADMOLEE SHARED ${SRC} )
endif()

if(LOOPLIB)
    target_link_libraries( HADMOLEE LOOPEXT )
endif()

target_link_libraries( HADMOLEE ${CUBATURELIB} )
target_link_libraries( HADMOLEE ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )

install(TARGETS HADMOLEE
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

if (BUILD_DICTIONARY)
    install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libHADMOLEE_rdict.pcm    
            DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

add_executable( hadMolee "cling/hadMolee.cpp")
target_link_libraries( hadMolee ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( hadMolee HADMOLEE )
target_link_libraries( hadMolee ${JSTYLELIB})
