#include "nonlinearfem/Kinematics2D.h" 
#include<iostream>
#include <cmath> 

int main() {
    using namespace NonlinearFEM ;

    constexpr double pi = 3.14159265358979323846;
    constexpr double theta = pi / 6.0;
    constexpr double tolerance = 1e-12;

    const double c = std::cos(theta);
    const double s = std::sin(theta);

    const Matrix2D displaceGradient {{
        {{c - 1.0 , -s}},
        {{s , c - 1.0}} 
    }} ;

    const Matrix2D deformationGradient =
        computeDeformationGradient(displaceGradient);

    const GreenLagrangeStrain Gstrain = 
        computeGreenLagrangeStrain (deformationGradient);
    
    std::cout << "Rigid rotation test:\n";

    std::cout << "Exx = " << Gstrain.exx << '\n';
    std::cout << "Eyy = " << Gstrain.eyy << '\n';
    std::cout << "Exy = " << Gstrain.exy << '\n';

    const bool passed =
        std::abs(Gstrain.exx) < tolerance &&
        std::abs(Gstrain.eyy) < tolerance &&
        std::abs(Gstrain.exy) < tolerance;

    if (!passed) {
        std::cerr << "Rigid rotation test failed.\n";
        return 1;
    }
        
    std::cout << "Rigid rotation test passed.\n";

    return 0;
}