set(CPP_SOURCES src/variant_matrix.cc src/fst.cc src/omega_max.cc src/polytable.cc src/summstats.cc src/windows_cpp.cc)

pybind11_add_module(_libsequence MODULE src/_libsequence.cc
    ${CPP_SOURCES})

target_link_libraries(_libsequence PRIVATE sequence)
