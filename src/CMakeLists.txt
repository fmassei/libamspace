if (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions(/W4)
elseif (CMAKE_COMPILER_IS_GNUCC)
    add_definitions("-Wall -Wextra")
endif()

#add_definitions("-g")
add_definitions("-O3 -mmmx -msse -msse2 -msse3")

add_definitions(-DBUILDING_LIBAMSPACE)
add_library (amspace SHARED angles.c coord.c date.c parseutils.c mtx.c
                            functions.c legendre.c interpolate.c integrate.c
                            linear.c leastSquare.c projections.c fk5.c
                            vsop87.c vsop87file.c elp82.c nutations.c
                            tle.c gravconst.c sat.c deep.c sgputils.c sgp4.c
                                sgp8.c)
