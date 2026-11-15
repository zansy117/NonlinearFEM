#include "nonlinearfem/Material.h" 
#include<iostream>


int main() {
    using namespace NonlinearFEM ;

    MaterialParameters material_;
    material_.poissonRatio = 0.3, material_.youngModulus = 200000;
    const MaterialParameters material = material_;

    GreenLagrangeStrain strain_;
    strain_.exx = 0.001, strain_.eyy = 0.002, strain_.exy = 0.0005;
    const GreenLagrangeStrain strain = strain_;

    SecondPKStress Stress =
        computeSecondPKStress(strain,material);
    
    std::cout<< Stress.sxx << " " << Stress.syy << " " << Stress.sxy << "\n";

    try {
        MaterialParameters invalidMaterial{
            200000.0,
            0.8
        };

        computeSecondPKStress(strain, invalidMaterial);

        std::cerr << "Error: invalid material was accepted.\n";
        return 1;
    }
    catch (const std::invalid_argument& error) {
        std::cout << "Invalid material test passed: "
                << error.what() << '\n';
    }
    
    return 0;
}