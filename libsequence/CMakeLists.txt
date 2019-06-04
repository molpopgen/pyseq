set(VARIANT_MATRIX_SOURCES src/variant_matrix.cc)

pybind11_add_module(_libsequence MODULE src/_libsequence.cc
    ${VARIANT_MATRIX_SOURCES})

target_link_libraries(_libsequence PRIVATE sequence)
