set(source_files
    constants.h
    file_names.cpp
    file_names.h
    utils.cpp
    utils.h
)

add_library(src ${source_files})

#Numerical library used
find_package(Armadillo REQUIRED)
#For parallelization
find_package(OpenMP REQUIRED)

#Couldve added OpenMP later, but whats the point
target_link_libraries(
    src
    PUBLIC
    armadillo
    OpenMP::OpenMP_CXX
)