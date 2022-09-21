# Build the executable to integrate with root cling
file(GLOB EXE_FILES "*.cpp")
foreach( exefile ${EXE_FILES} )
    get_filename_component( exename ${exefile} NAME_WE)
    add_executable( ${exename} ${exefile} )
    target_link_libraries( ${exename} ROOT::Core ROOT::MathMore ROOT::RIO ROOT::Rint )
    target_link_libraries( ${exename} EE3BODY )
    target_link_libraries( ${exename} ${JSTYLELIB})
    
    set_property(TARGET ${exename} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries( ${exename} ${LOOPLIB} )

endforeach( exefile ${EXE_FILES} )