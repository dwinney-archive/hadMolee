# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "core/*.cpp" "analysis/*.cpp" "cubature/*.c")


foreach(file_path ${SRC})
    message(${file_path})
endforeach()

include_directories(core)
include_directories(cubature)
include_directories(analysis)

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


if(LOOPLIB)
    target_link_libraries( HADMOLEE LOOPEXT )
endif()

target_link_libraries( HADMOLEE ROOT::Core ROOT::MathMore )

install(TARGETS HADMOLEE
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )


install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libHADMOLEE_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable( hadMolee "cling/hadMolee.cpp")
target_link_libraries( hadMolee ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
target_link_libraries( hadMolee HADMOLEE )
