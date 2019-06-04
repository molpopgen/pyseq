set(CPP_SOURCES src/variant_matrix.cc src/fst.cc src/omega_max.cc src/polytable.cc src/summstats.cc src/windows_cpp.cc)
file(GLOB LIBSEQ_SOURCES src/libsequence/src/*.cc src/libsequence/src/Seq/*.cc
    src/libsequence/src/variant_matrix/*.cc 
    src/libsequence/src/summstats/*.cc 
    src/libsequence/src/summstats_deprecated/*.cc)
get_filename_component(full_path_test_cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/libsequence/src/libsequenceConfig.cc ABSOLUTE)
list(REMOVE_ITEM LIBSEQ_SOURCES "${full_path_test_cpp}")
get_filename_component(full_path_test_cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/libsequence/src/PolyTableManip.cc ABSOLUTE)
list(REMOVE_ITEM LIBSEQ_SOURCES "${full_path_test_cpp}")
pybind11_add_module(_libsequence MODULE src/_libsequence.cc
    ${CPP_SOURCES}
    ${LIBSEQ_SOURCES})

# target_link_libraries(_libsequence PRIVATE sequence)
