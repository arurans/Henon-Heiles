set(
    problem_files
    compute.cpp
    compute.h
    hamiltonian.cpp
    hamiltonian.h
    poincare.cpp
    poincare.h
)

add_library(problems ${problem_files})

target_link_libraries(
    problems
    methods
)