cmake_minimum_required(VERSION 3.16.3)

project(hhp LANGUAGES CXX) # Hénon Heiles (parallelized)

# -O1 seems to be faster than both -O2 and -O3 in this case
# and my system
# CPU: Ryzen 7 5600x (not overclocked at the time of computing)
# RAM: DDR4 32GB 3600 Mhz
set(CMAKE_CXX_FLAGS -O1)


add_subdirectory(src/problems)
add_subdirectory(src/methods)
add_subdirectory(src)


add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(
    ${PROJECT_NAME}
    problems
)

# Can uncomment the lines below to run the file
# automatically after building it
# It works fine on my system, but I have no idea
# if doing such things might cause some unexpected bugs,
# given a small change in the code or another system...
# It's C++ after all (^:

# add_custom_command(
#     TARGET 
#         ${PROJECT_NAME}
#     POST_BUILD
#     COMMENT
#         "Running main.cpp"
#     COMMAND 
#         ${PROJECT_NAME}
# )

#get_filename_component(par ../ ABSOLUTE)