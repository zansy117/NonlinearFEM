#include "nonlinearfem/Assembly.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using namespace NonlinearFEM;
namespace {
void near(double actual,double expected,const char* label) {
    if (!std::isfinite(actual) || std::abs(actual-expected)>1e-11*std::max(1.0,std::abs(expected)))
        throw std::runtime_error(label);
}
template<class Exception,class Function>
void rejects(Function function,const char* label) {
    try { function(); } catch(const Exception&) { return; }
    throw std::runtime_error(label);
}
Q4Mesh pairMesh() {
    return {{{0,0},{1,0},{2,0},{0,1},{1,1},{2,1}},{{0,1,4,3},{1,2,5,4}}};
}
}
int main() {
    try {
        const Q4NodalValues square={{{0,0},{1,0},{1,1},{0,1}}};
        const auto local=computeQ4LumpedMass(square,1000,0.2);
        for(double mass:local) near(mass,50,"Single-element nodal mass");
        const Q4Mesh single={{{0,0},{1,0},{1,1},{0,1}},{{0,1,2,3}}};
        const auto singleGlobal=assembleGlobalLumpedMass(single,1000,0.2);
        for(std::size_t i=0;i<8;++i) near(singleGlobal.at(i),local[i],"Single-element scatter");
        const auto mesh=pairMesh();
        const auto mass=assembleGlobalLumpedMass(mesh,1000,0.2);
        if(mass.size()!=12) throw std::runtime_error("Mass vector size");
        double sumX=0,sumY=0;
        for(std::size_t n=0;n<6;++n) {
            const double expected=(n==1 || n==4)?100:50;
            near(mass[2*n],expected,"Shared-node X mass");
            near(mass[2*n+1],expected,"Shared-node Y mass");
            sumX+=mass[2*n];sumY+=mass[2*n+1];
        }
        near(sumX,400,"Total X mass");near(sumY,400,"Total Y mass");
        const auto twiceDensity=assembleGlobalLumpedMass(mesh,2000,0.2);
        const auto twiceThickness=assembleGlobalLumpedMass(mesh,1000,0.4);
        for(std::size_t i=0;i<mass.size();++i) {
            near(twiceDensity[i],2*mass[i],"Density scaling");
            near(twiceThickness[i],2*mass[i],"Thickness scaling");
        }
        rejects<std::invalid_argument>([]{(void)assembleGlobalLumpedMass(Q4Mesh{},1000,0.2);},"Empty mesh accepted");
        auto bad=mesh;bad.elements.clear();
        rejects<std::invalid_argument>([&]{(void)assembleGlobalLumpedMass(bad,1000,0.2);},"Empty connectivity accepted");
        bad=mesh;bad.elements[0][0]=99;
        rejects<std::out_of_range>([&]{(void)assembleGlobalLumpedMass(bad,1000,0.2);},"Invalid index accepted");
        bad=mesh;bad.elements[0][1]=bad.elements[0][0];
        rejects<std::invalid_argument>([&]{(void)assembleGlobalLumpedMass(bad,1000,0.2);},"Duplicate node accepted");
        bad=mesh;bad.nodes[0][0]=std::numeric_limits<double>::quiet_NaN();
        rejects<std::invalid_argument>([&]{(void)assembleGlobalLumpedMass(bad,1000,0.2);},"Non-finite coordinate accepted");
        for(double thickness:{0.0,-1.0,std::numeric_limits<double>::infinity()})
            rejects<std::invalid_argument>([&]{(void)assembleGlobalLumpedMass(mesh,1000,thickness);},"Invalid thickness accepted");
        // A distorted multi-element mesh and a non-affine displacement field.
        auto distorted=mesh;
        distorted.nodes[4]={1.1,1.2};distorted.nodes[5]={2.2,1.0};
        const GlobalDisplacement displacement={0.01,-0.02,0.03,0.01,0.06,-0.01,-0.01,0.04,0.05,0.07,0.08,0.02};
        MaterialParameters material{};material.youngModulus=210000;material.poissonRatio=0.3;
        const auto force=assembleInternalForce(distorted,displacement,material,0.2);
        const auto system=assembleGlobalSystem(distorted,displacement,material,0.2);
        if(force.size()!=system.internalForce.size()) throw std::runtime_error("Internal force size");
        for(std::size_t i=0;i<force.size();++i)
            near(force[i],system.internalForce[i],"Force-only assembly mismatch");
        std::cout<<"PASS: element mass, shared-node scatter, total mass, scaling, geometry rejection and force-only assembly.\n";
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
    return 0;
}