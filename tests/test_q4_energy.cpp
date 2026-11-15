#include "nonlinearfem/Assembly.h"
#include "nonlinearfem/ExplicitDynamics.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace NonlinearFEM;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
const Q4NodalValues square={{{0,0},{1,0},{1,1},{0,1}}};
MaterialParameters material(){MaterialParameters m{};m.youngModulus=1000;m.poissonRatio=.3;return m;}
void gradientTest(){
    const Q4NodalValues X={{{0,0},{1.2,.1},{1.1,1.2},{-.1,.9}}};
    const Q4NodalValues u={{{.01,-.02},{.04,.01},{.06,.07},{-.02,.03}}};
    const auto m=material();
    const auto force=computeQ4InternalForce(X,u,m,.2);
    const double h=1e-6;
    double error=0,reference=0;
    for(std::size_t i=0;i<8;++i){
        auto plus=u,minus=u;plus[i/2][i%2]+=h;minus[i/2][i%2]-=h;
        const double derivative=(computeQ4StrainEnergy(X,plus,m,.2)-computeQ4StrainEnergy(X,minus,m,.2))/(2*h);
        error+=(derivative-force[i])*(derivative-force[i]);reference+=force[i]*force[i];
    }
    const double relative=std::sqrt(error/reference);
    std::cout<<"Distorted Q4 energy-gradient relative error="<<relative<<'\n';
    require(std::isfinite(relative)&&relative<1e-7,"Energy derivative does not match internal force");
    require(computeQ4StrainEnergy(square,Q4NodalValues{},m,.2)==0,"Zero-state energy");
    // Analytic homogeneous simple shear, including finite-strain normal energy.
    const double gamma=.1;
    Q4NodalValues shear{};for(std::size_t a=0;a<4;++a)shear[a][0]=gamma*square[a][1];
    const double eyy=.5*gamma*gamma;
    const double expected=.2*.5*(m.youngModulus/(1-m.poissonRatio*m.poissonRatio)*eyy*eyy
        +m.youngModulus/(2*(1+m.poissonRatio))*gamma*gamma);
    const double energy=computeQ4StrainEnergy(square,shear,m,.2);
    require(std::isfinite(energy)&&std::abs(energy-expected)<1e-12,"Analytic finite-shear energy mismatch");
    require(std::abs(computeQ4StrainEnergy(square,shear,m,.4)-2*energy)<1e-12,"Energy thickness scaling");
}
double energyHistory(int refinement){
    const Q4Mesh mesh={{{0,0},{1,0},{1,1},{0,1}},{{0,1,2,3}}};
    const auto m=material();const auto mass=assembleGlobalLumpedMass(mesh,1,.2);
    const std::vector<bool> fixed={true,true,false,false,false,false,true,true};
    ExplicitState s{};s.displacement.assign(8,0);s.velocity.assign(8,0);
    s.displacement[2]=s.displacement[4]=1e-4;
    const InternalForceFunction internal=[&](const DynamicVector& u){return assembleInternalForce(mesh,u,m,.2);};
    const ExternalForceFunction external=[](double){return DynamicVector(8,0);};
    const auto totalEnergy=[&](const ExplicitState& state){
        Q4NodalValues u{};double kinetic=0;
        for(std::size_t i=0;i<8;++i){u[i/2][i%2]=state.displacement[i];kinetic+=.5*mass[i]*state.velocity[i]*state.velocity[i];}
        return kinetic+computeQ4StrainEnergy(square,u,m,.2);
    };
    InitializeExplicitAcceleration(s,mass,fixed,internal,external);
    const double initial=totalEnergy(s);require(initial>0,"Positive initial energy required");
    double maximum=0;
    for(int i=0;i<2000*refinement;++i){
        advanceCentralDifference(s,mass,fixed,1e-4/refinement,internal,external);
        const double energy=totalEnergy(s);
        require(std::isfinite(energy)&&energy>=0,"Invalid total energy");
        maximum=std::max(maximum,std::abs(energy-initial)/initial);
    }
    std::cout<<"dt="<<1e-4/refinement<<", maximum relative total-energy error="<<maximum<<'\n';
    require(maximum<1e-4,"Excessive free-vibration energy error");
    return maximum;
}
}
int main(){try{
    gradientTest();
    const double a=energyHistory(1),b=energyHistory(2),c=energyHistory(4);
    require(a>0&&b>0&&c>0&&a/b>3.5&&a/b<4.5&&b/c>3.5&&b/c<4.5,"Energy error does not decrease quadratically");
    std::cout<<"PASS: analytic strain energy, energy-force consistency and free-vibration energy refinement.\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}return 0;}