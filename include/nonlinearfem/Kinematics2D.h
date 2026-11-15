#pragma once

#include <array>
#include <cstddef>

namespace NonlinearFEM {

    constexpr std::size_t kSpatialDimension = 2;

    using Matrix2D =
        std::array<std::array<double, 2>, 2 >;
    

    struct GreenLagrangeStrain {
        double exx = 0.0 ;
        double eyy = 0.0 ;
        double exy = 0.0 ;
    } ;
     
    [[nodiscard]]
    Matrix2D computeDeformationGradient(
        const Matrix2D& displacementGradient
    );

    [[nodiscard]]
    GreenLagrangeStrain computeGreenLagrangeStrain(
        const Matrix2D& deformationGradient
    ); 


} // namespace NonlinearFEM