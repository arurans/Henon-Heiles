set(source_files
    constants.h
    storage_info.cpp
    storage_info.h
    utils.cpp
    utils.h
)

add_library(src ${source_files})

#Numerical library used
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#For parallelization
find_package(OpenMP REQUIRED)

#Couldve added OpenMP later, but whats the point
target_link_libraries(
    src
    PUBLIC
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
)