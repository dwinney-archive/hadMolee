# BUILD LIBRARY FROM LOCAL FiLES
file(GLOB SRC "*.cpp")

ROOT_GENERATE_DICTIONARY(G__EE3BODY
        reaction_kinematics.hpp
        amplitude.hpp
        propagator.hpp
        LINKDEF EE3BODYLinkDef.h
        )

add_library( EE3BODY SHARED ${SRC} G__EE3BODY.cxx )
target_link_libraries( EE3BODY ROOT::Core ROOT::MathMore )
target_link_libraries( EE3BODY Boost::boost)

install(TARGETS EE3BODY
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libEE3BODY_rdict.pcm    
         DESTINATION ${CMAKE_INSTALL_LIBDIR})