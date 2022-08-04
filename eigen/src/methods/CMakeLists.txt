set(
    method_files
    kahans.cpp
    kahans.h
    rk4.cpp
    rk4.h
    sb.cpp
    sb.h
    sv.cpp
    sv.h
)

add_library(methods ${method_files})

target_link_libraries(
    methods
    src
)