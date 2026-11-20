#pragma once

#include <array>
#include <cstddef>

namespace NonlinearFEM {

    constexpr std::size_t kSpatialDimension_3D = 3;

    using Matrix3D = std::array<std::array<double, 3>, 3>;
    

    struct GreenLagrangeStrain_3D {
        double exx = 0.0;
        double eyy = 0.0;
        double ezz = 0.0;
        double exy = 0.0;
        double exz = 0.0;
        double eyz = 0.0;
    } ;
     
    [[nodiscard]]
    Matrix3D computeDeformationGradient_3D(
        const Matrix3D& displacementGradient
    );

    [[nodiscard]]
    GreenLagrangeStrain_3D computeGreenLagrangeStrain_3D(
        const Matrix3D& deformationGradient
    ); 


} // namespace NonlinearFEM