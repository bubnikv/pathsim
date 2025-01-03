cmake_minimum_required(VERSION 3.8)
project(PathSim)

include_directories("${PROJECT_SOURCE_DIR}")

set(PathSimSources
    cmplx.h
    Delay.cpp
    Delay.h
    FilterTables.h
    GaussFIR.cpp
    GaussFIR.h
    main.cpp
    NoiseGen.cpp
    NoiseGen.h
    Path.cpp
    Path.h
    PathSimParams.cpp
    PathSimParams.h
    PathSimProcessor.cpp
    PathSimProcessor.h
    )

add_executable(pathsim ${PathSimSources})

#install(TARGETS pathsim RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
