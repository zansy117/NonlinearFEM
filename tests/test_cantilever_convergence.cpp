#include "nonlinearfem/LoadStepping.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace NonlinearFEM;
namespace {
constexpr double length = 4.0, height = 1.0, thickness = 0.2;
constexpr double young = 1.0e5, poisson = 0.3;
struct Response { double ux, uy; };

Response solve(int nx, int ny, double force) {
    Q4Mesh mesh{};
    for (int j=0; j<=ny; ++j)
        for (int i=0; i<=nx; ++i)
            mesh.nodes.push_back({length*i/nx, height*j/ny});
    for (int j=0; j<ny; ++j)
        for (int i=0; i<nx; ++i) {
            const auto a=static_cast<std::size_t>(j*(nx+1)+i);
            const auto b=a+static_cast<std::size_t>(nx+1);
            mesh.elements.push_back({a,a+1,b+1,b});
        }
    DisplacementConstraints constraints;
    NodalLoads loads;
    for (int j=0; j<=ny; ++j) {
        const auto left=static_cast<std::size_t>(j*(nx+1));
        constraints.push_back({left,0,0.0});
        constraints.push_back({left,1,0.0});
        // Consistent nodal forces for uniform traction on the reference edge.
        const double weight=(j==0 || j==ny) ? 0.5 : 1.0;
        loads.push_back({left+static_cast<std::size_t>(nx),1,-force*weight/ny});
    }
    MaterialParameters material{};
    material.youngModulus=young;
    material.poissonRatio=poisson;
    LoadSteppingOptions options{};
    options.initialIncrement=0.25;
    options.maximumIncrement=0.25;
    options.newtonOptions.maxIterations=40;
    options.newtonOptions.residualTolerance=1.0e-9;
    const auto result=solveLoadStepping(mesh,material,loads,constraints,
        GlobalDisplacement(2*mesh.nodes.size(),0.0),options,thickness);
    if (!result.converged || std::abs(result.finalLoadFactor-1.0)>1.0e-12)
        throw std::runtime_error("Load stepping failed: "+result.failureReason);
    const auto assembly=assembleGlobalSystem(mesh,result.displacement,material,thickness);
    const auto external=assembleExternalForce(mesh.nodes.size(),loads);
    double residualSquared=0.0, rx=0.0, ry=0.0, constraintError=0.0;
    for (std::size_t node=0; node<mesh.nodes.size(); ++node) {
        for (std::size_t d=0; d<2; ++d) {
            const auto k=2*node+d;
            if (!std::isfinite(result.displacement[k]))
                throw std::runtime_error("Non-finite displacement");
            const double residual=assembly.internalForce[k]-external[k];
            if (node%static_cast<std::size_t>(nx+1)==0) {
                constraintError=std::max(constraintError,std::abs(result.displacement[k]));
                if (d==0) rx+=residual; else ry+=residual;
            } else residualSquared+=residual*residual;
        }
    }
    const double relativeResidual=std::sqrt(residualSquared)/force;
    const double balance=std::hypot(rx,ry-force)/force;
    if (!(relativeResidual<1.0e-6 && balance<1.0e-6 && constraintError<1.0e-10))
        throw std::runtime_error("Equilibrium or constraint verification failed");
    // Interpolate the reference-edge midpoint for odd ny as well as even ny.
    const auto low=static_cast<std::size_t>((ny/2)*(nx+1)+nx);
    const auto high=static_cast<std::size_t>(((ny+1)/2)*(nx+1)+nx);
    const Response response={0.5*(result.displacement[2*low]+result.displacement[2*high]),
                            0.5*(result.displacement[2*low+1]+result.displacement[2*high+1])};
    std::cout<<nx<<','<<ny<<','<<force<<','<<response.ux<<','<<response.uy
             <<','<<relativeResidual<<','<<balance<<'\n';
    if (!(response.uy<0.0)) throw std::runtime_error("Wrong deflection direction");
    return response;
}
}
int main() {
    try {
        const std::array<int,4> nx={4,8,16,24}, ny={1,2,4,6};
        const std::array<double,5> forces={0.1,12.5,25.0,37.5,50.0};
        std::array<std::array<Response,5>,4> responses{};
        std::cout<<std::scientific<<std::setprecision(10)
                 <<"nx,ny,downward_force,tip_ux,tip_uy,free_residual,force_balance\n";
        for (std::size_t m=0; m<nx.size(); ++m) {
            double previous=0.0;
            for (std::size_t p=0; p<forces.size(); ++p) {
                responses[m][p]=solve(nx[m],ny[m],forces[p]);
                if (!(-responses[m][p].uy>previous))
                    throw std::runtime_error("Non-monotone load-deflection curve");
                previous=-responses[m][p].uy;
            }
        }
        // Beam theory is only an approximate small-load reference for this 2D solid.
        const double inertia=thickness*height*height*height/12.0;
        const double bending=forces[0]*length*length*length/(3.0*young*inertia);
        const double shearModulus=young/(2.0*(1.0+poisson));
        const double shear=forces[0]*length/((5.0/6.0)*shearModulus*thickness*height);
        const double beam=bending+shear;
        const double beamError=std::abs(-responses[3][0].uy-beam)/beam;
        std::cout<<"Small-load Timoshenko estimate: "<<beam
                 <<"; finest-mesh relative difference: "<<beamError<<'\n';
        if (!(beamError<0.10)) throw std::runtime_error("Small-load beam comparison failed");
        for (std::size_t p=0; p<forces.size(); ++p) {
            const double coarseChange=std::abs(responses[2][p].uy-responses[1][p].uy)
                                      /std::abs(responses[2][p].uy);
            const double fineChange=std::abs(responses[3][p].uy-responses[2][p].uy)
                                    /std::abs(responses[3][p].uy);
            std::cout<<"Force "<<forces[p]<<": relative mesh changes "
                     <<coarseChange<<" -> "<<fineChange<<'\n';
            if (!(fineChange<coarseChange && fineChange<0.05))
                throw std::runtime_error("Mesh refinement criterion failed");
        }
        std::cout<<"Cantilever mesh refinement verification passed.\n";
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n'; return 1;
    }
    return 0;
}