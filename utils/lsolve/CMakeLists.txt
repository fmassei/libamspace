if (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions(/W4)
elseif (CMAKE_COMPILER_IS_GNUCC)
    add_definitions("-Wall -Wextra")
endif()

include_directories("${PROJECT_SOURCE_DIR}/src")
link_directories("${PROJECT_SOURCE_DIR}/src")

#add_definitions("-g")

add_executable(lsolve main.c)
target_link_libraries(lsolve amspace m)
