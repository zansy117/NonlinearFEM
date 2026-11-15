
cmake_minimum_required(VERSION 3.16)

project(NonlinearFEM LANGUAGES CXX)

include(CTest)

# ==================================
# C++ standard
# ==================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ==================================
# Main executable
# ==================================

add_executable(NonlinearFEM
    src/main.cpp
    src/ModelIO.cpp
    src/Solver.cpp
    src/BoundaryConditions.cpp
    src/LoadStepping.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
    src/ElementQ4.cpp
    src/Assembly.cpp
)

# ==================================
# Kinematics test
# ==================================

add_executable(test_kinematics
    tests/test_kinematics.cpp
    src/Kinematics2D.cpp
)

# ==================================
# Green-Lagrange strain test
# ==================================

add_executable(test_Gstrain
    tests/test_Gstrain.cpp
    src/Kinematics2D.cpp
)

# ==================================
# Material test
# ==================================

add_executable(test_material
    tests/test_material.cpp
    src/Material.cpp
)

# ==================================
# Q4 element test
# ==================================

add_executable(test_q4
    tests/test_q4.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

# ==================================
# Include directories
# ==================================

add_executable(test_assembly
    tests/test_assembly.cpp
    src/Assembly.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_boundary
    tests/test_boundary.cpp
    src/BoundaryConditions.cpp
    src/Assembly.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_solver
    tests/test_solver.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_solver_2x2
    tests/test_solver_2x2.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_tangent_stiffness
    tests/test_tangent_stiffness.cpp
    src/Assembly.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_load_stepping
    tests/test_load_stepping.cpp
    src/LoadStepping.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_load_stepping_cutback
    tests/test_load_stepping_cutback.cpp
    src/LoadStepping.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_load_stepping_displacement
    tests/test_load_stepping_displacement.cpp
    src/LoadStepping.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_q4_patch_multi
    tests/test_q4_patch_multi.cpp
    src/LoadStepping.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_q4_patch_stress
    tests/test_q4_patch_stress.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_q4_patch_distorted
    tests/test_q4_patch_distorted.cpp
    src/LoadStepping.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_tangent_stiffness_distorted
    tests/test_tangent_stiffness_distorted.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_global_assembly_verification
    tests/test_global_assembly_verification.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
    src/ElementQ4.cpp
    src/Assembly.cpp
)

add_executable(test_uniaxial_tension
    tests/test_uniaxial_tension.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_cantilever_convergence
    tests/test_cantilever_convergence.cpp
    src/LoadStepping.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_model_io
    tests/test_model_io.cpp
    src/ModelIO.cpp
    src/LoadStepping.cpp
    src/Solver.cpp
    src/Assembly.cpp
    src/BoundaryConditions.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_lumped_mass
    tests/test_lumped_mass.cpp
    src/Assembly.cpp
    src/ElementQ4.cpp
    src/Kinematics2D.cpp
    src/Material.cpp
)

add_executable(test_explicit_dynamics
    tests/test_explicit_dynamics.cpp
    src/ExplicitDynamics.cpp
)

set(NONLINEARFEM_TARGETS
    NonlinearFEM
    test_explicit_dynamics
    test_lumped_mass
    test_model_io
    test_cantilever_convergence
    test_uniaxial_tension
    test_kinematics
    test_Gstrain
    test_material
    test_q4
    test_assembly
    test_boundary
    test_solver
    test_solver_2x2
    test_tangent_stiffness
    test_load_stepping
    test_load_stepping_cutback
    test_load_stepping_displacement
    test_q4_patch_multi
    test_q4_patch_stress
    test_q4_patch_distorted
    test_tangent_stiffness_distorted
    test_global_assembly_verification
)

foreach(target IN LISTS NONLINEARFEM_TARGETS)

    target_include_directories(${target}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

endforeach()

# ==================================
# Compiler warnings
# ==================================

if(MSVC)

    foreach(target IN LISTS NONLINEARFEM_TARGETS)
        target_compile_options(${target} PRIVATE /W4)
    endforeach()

else()

    foreach(target IN LISTS NONLINEARFEM_TARGETS)
        target_compile_options(${target} PRIVATE
            -Wall
            -Wextra
            -Wpedantic
        )
    endforeach()

endif()

# Register the existing test executables with CTest.
if(BUILD_TESTING)
    foreach(target IN LISTS NONLINEARFEM_TARGETS)
        if(target MATCHES "^test_")
            if(target STREQUAL "test_model_io")
                add_test(NAME ${target} COMMAND ${target} "${CMAKE_CURRENT_SOURCE_DIR}/examples/tension.nfem")
            else()
                add_test(NAME ${target} COMMAND ${target})
            endif()
        endif()
    endforeach()
endif()

# The local slider interface uses only Python's standard library.
if(BUILD_TESTING)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME test_gui_model
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gui_model.py" "$<TARGET_FILE:NonlinearFEM>")
    endif()
endif()