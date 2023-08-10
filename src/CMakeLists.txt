# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "core/*.cpp" "analysis/*.cpp" "cubature/*.c")

include_directories(core)
include_directories(cubature)
include_directories(analysis)
include_directories(${LOOPTOOLS_INC})

# Build LoopTools
add_library(LOOPTOOLS SHARED ${STATIC_LOOPTOOLS})
set_target_properties(LOOPTOOLS PROPERTIES LINKER_LANGUAGE Fortran )
target_link_libraries(LOOPTOOLS PUBLIC
        -Wl,--whole-archive
        ${STATIC_LOOPTOOLS}
        -Wl,--no-whole-archive
        )

# Build hadMolee and ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__HADMOLEE
        kinematics.hpp
        amplitude.hpp
        lineshape.hpp
        breit_wigner.hpp
        triangle.hpp
        box.hpp
        molecule.hpp
        fitter.hpp
        plotter.hpp
        plot.hpp
        LINKDEF HADMOLEELinkDef.h
        )
add_library( HADMOLEE SHARED ${SRC} G__HADMOLEE.cxx )
target_link_libraries( HADMOLEE ROOT::Core ROOT::MathMore )

# INSTALLATION 
install(TARGETS LOOPTOOLS
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

install(TARGETS HADMOLEE
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libHADMOLEE_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable( hadMolee "cling/hadMolee.cpp")
target_link_libraries( hadMolee ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( hadMolee HADMOLEE )
