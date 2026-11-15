#include "nonlinearfem/Kinematics2D.h" 
#include<iostream>
#include <cmath> 

int main() {
    using namespace NonlinearFEM ;

    const Matrix2D displacementGradient {{
        {{0.1, 0.2}},
        {{0.3, 0.4}} 
    }} ;

    const Matrix2D expected{{
        {{1.1, 0.2}},
        {{0.3, 1.4}}
    }};

    const Matrix2D deformationGradient = 
        computeDeformationGradient(displacementGradient) ;

    constexpr double tolerance = 1e-12 ;

    bool passed = true ;

    for(std::size_t i = 0; i < kSpatialDimension ; ++i ) {
        for(std::size_t j = 0; j < kSpatialDimension ; ++j ) {
            const double error = std::abs(
                deformationGradient[i][j] - expected[i][j]
            );

            if (error > tolerance) {
                passed = false;
            }

            std::cout << deformationGradient[i][j] << " ";
        }
        std::cout <<"\n";
    }    

    if (passed) {
        std::cout << "Test passed.\n";
        return 0;
    }

    std::cerr << "Test failed.\n";
    return 1;
    
}
