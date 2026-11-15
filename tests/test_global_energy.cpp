#include "nonlinearfem/Assembly.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace NonlinearFEM;
namespace {
void near(double a,double b,double tolerance,const char* message){if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tolerance*std::max(1.0,std::abs(b)))throw std::runtime_error(message);}
}
int main(){try{
    Q4Mesh mesh={{{0,0},{1,0},{2.1,.1},{0,1},{1.1,1.2},{2.2,1.1}},{{0,1,4,3},{1,2,5,4}}};
    const GlobalDisplacement u={.01,-.02,.03,.01,.07,-.01,-.02,.04,.06,.08,.09,.03};
    MaterialParameters m{};m.youngModulus=210000;m.poissonRatio=.3;
    const double thickness=.2;
    const double energy=assembleGlobalEnergy(mesh,u,m,thickness);
    double sum=0;
    for(const auto& e:mesh.elements){
        Q4NodalValues X{},local{};
        for(std::size_t a=0;a<4;++a){X[a]=mesh.nodes[e[a]];local[a]={u[2*e[a]],u[2*e[a]+1]};}
        sum+=computeQ4StrainEnergy(X,local,m,thickness);
    }
    near(energy,sum,1e-12,"Element energy sum mismatch");
    auto reversed=mesh;std::reverse(reversed.elements.begin(),reversed.elements.end());
    near(assembleGlobalEnergy(reversed,u,m,thickness),energy,1e-12,"Element ordering changed energy");
    near(assembleGlobalEnergy(mesh,GlobalDisplacement(12,0),m,thickness),0,1e-14,"Zero-state energy");
    near(assembleGlobalEnergy(mesh,u,m,2*thickness),2*energy,1e-12,"Thickness scaling");
    const auto force=assembleInternalForce(mesh,u,m,thickness);
    double error=0,norm=0;
    const double h=1e-6;
    for(std::size_t i=0;i<u.size();++i){
        auto plus=u,minus=u;plus[i]+=h;minus[i]-=h;
        const double derivative=(assembleGlobalEnergy(mesh,plus,m,thickness)-assembleGlobalEnergy(mesh,minus,m,thickness))/(2*h);
        error+=(derivative-force[i])*(derivative-force[i]);norm+=force[i]*force[i];
    }
    const double relative=std::sqrt(error/norm);
    std::cout<<"Multi-element energy="<<energy<<", gradient relative error="<<relative<<'\n';
    if(!std::isfinite(relative)||relative>1e-7)throw std::runtime_error("Global energy derivative disagrees with internal force");
    // Independent homogeneous finite-shear energy on a rectangular two-element mesh.
    mesh.nodes={{0,0},{1,0},{2,0},{0,1},{1,1},{2,1}};
    GlobalDisplacement shear(12,0);const double gamma=.1;
    for(std::size_t n=0;n<6;++n)shear[2*n]=gamma*mesh.nodes[n][1];
    const double eyy=.5*gamma*gamma;
    const double density=.5*(m.youngModulus/(1-m.poissonRatio*m.poissonRatio)*eyy*eyy
        +m.youngModulus/(2*(1+m.poissonRatio))*gamma*gamma);
    near(assembleGlobalEnergy(mesh,shear,m,thickness),density*2*thickness,1e-12,"Analytic multi-element shear energy");
    std::cout<<"PASS: element sum, ordering, zero energy, thickness scaling, global force gradient and analytic shear.\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}return 0;}