if (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions(/W4)
elseif (CMAKE_COMPILER_IS_GNUCC)
    add_definitions("-Wall -Wextra")
endif()

include_directories("${PROJECT_SOURCE_DIR}/src")
link_directories("${PROJECT_SOURCE_DIR}/src")

add_definitions("-O3 -mmmx -msse -msse2 -msse3")

add_executable(sgpperformace main.c)
target_link_libraries(sgpperformace amspace m pthread)
