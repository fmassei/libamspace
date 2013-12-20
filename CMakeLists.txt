cmake_minimum_required (VERSION 2.6)
project (libamspace)

set (LIBAMSPACE_VERSION_MAJOR 0)
set (LIBAMSPACE_VERSION_MINOR 1)
set (LIBAMSPACE_VERSION_REV \"1beta\")

configure_file (
    "${PROJECT_SOURCE_DIR}/src/config.h.in"
    "${PROJECT_BINARY_DIR}/src/config.h"
    )

include_directories("${PROJECT_BINARY_DIR}/src")

add_subdirectory (src)
add_subdirectory (utils)
add_subdirectory (test)

# testing
include (CTest)

add_test (sgprun utils/sgp/sgp -v)

macro (test_sgp_pre tname sname ver)
    file (WRITE ${CMAKE_BINARY_DIR}/tmptst${tname}.sh
        "cat ${PROJECT_SOURCE_DIR}/test/data/${sname}.tle | \\
            utils/sgp/sgp -m ${ver} | \\
            test/sgpverify/sgpverify - \\
                ${PROJECT_SOURCE_DIR}/test/data/${sname}.dat.orig.sgp${ver} s")
    add_test(${tname} sh ${CMAKE_BINARY_DIR}/tmptst${tname}.sh)
    set_tests_properties (${tname}
        PROPERTIES PASS_REGULAR_EXPRESSION
        "^-?0[.][0-1][^ ]* -?0[.]000[0-9]*")
endmacro (test_sgp_pre)

test_sgp_pre ("SGP4" "88888" "4")
test_sgp_pre ("SDP4" "11801" "4")
test_sgp_pre ("SGP8" "88888" "8")
test_sgp_pre ("SDP8" "11801" "8")

add_test (VSOP87 test/vsop87test/vsop87test ${PROJECT_SOURCE_DIR}/data/vsop87
            ${PROJECT_SOURCE_DIR}/test/data/vsop87.chk)
add_test (ELP82 utils/elp82/elp82 ${PROJECT_SOURCE_DIR}/data/elp82b/ test)

add_test (mathfunctions test/testmath/testmath f)
add_test (mathlegendre test/testmath/testmath l)
add_test (mathinterpolate test/testmath/testmath i)
add_test (mathintegrate test/testmath/testmath n)
add_test (mathlinear test/testmath/testmath e)
add_test (mathleastsquare test/testmath/testmath s)

add_test (projections test/testprojections/testprojections)

file (WRITE ${CMAKE_BINARY_DIR}/sgpcatalogtest.sh
    "cat ${PROJECT_SOURCE_DIR}/test/data/SGP4-VER.TLE | \\
        utils/sgp/sgp -t | \\
        test/sgpverify/sgpverify - \\
        ${PROJECT_SOURCE_DIR}/test/data/SGP4-VER.dat.orig")

file (WRITE ${CMAKE_BINARY_DIR}/sgp8catalogtest.sh
    "cat ${PROJECT_SOURCE_DIR}/test/data/SGP4-VER.TLE | \\
        utils/sgp/sgp -t -m 8 | \\
        test/sgpverify/sgpverify - \\
        ${PROJECT_SOURCE_DIR}/test/data/SGP4-VER.dat.orig")
