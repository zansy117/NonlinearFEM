#include "nonlinearfem/ExplicitDynamics.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace NonlinearFEM;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void near(double a,double b,double tolerance,const char* message){require(std::isfinite(a)&&std::abs(a-b)<=tolerance,message);}
const InternalForceFunction spring=[](const DynamicVector& u){return DynamicVector{100*u[0]};};
const ExternalForceFunction zero=[](double){return DynamicVector{0};};
ExplicitState initial(){ExplicitState s{};s.displacement={.01};s.velocity={0};return s;}
void unchanged(const ExplicitState& a,const ExplicitState& b){require(a.time==b.time&&a.displacement==b.displacement&&a.velocity==b.velocity&&a.acceleration==b.acceleration,"Failed step changed state");}
template<class F> void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,"Invalid input accepted");}
double oscillator(int steps){
    auto s=initial();const double dt=1.0/steps;
    InitializeExplicitAcceleration(s,{1},{false},spring,zero);
    near(s.acceleration.at(0),-1,1e-14,"Initial acceleration");
    double error=0,maxEnergyError=0;
    for(int i=0;i<steps;++i){
        advanceCentralDifference(s,{1},{false},dt,spring,zero);
        error=std::max(error,std::abs(s.displacement[0]-.01*std::cos(10*s.time)));
        const double energy=.5*s.velocity[0]*s.velocity[0]+50*s.displacement[0]*s.displacement[0];
        maxEnergyError=std::max(maxEnergyError,std::abs(energy-.005)/.005);
        near(s.acceleration[0],-100*s.displacement[0],1e-12,"Acceleration equilibrium");
    }
    near(s.time,1,1e-12,"Final time");
    require(error<5e-5,"Oscillator displacement error");
    require(maxEnergyError<.003,"Oscillator energy error");
    near(s.velocity[0],-.1*std::sin(10*s.time),5e-4,"Oscillator velocity error");
    std::cout<<"dt="<<dt<<", maximum displacement error="<<error<<", relative energy error="<<maxEnergyError<<'\n';
    return error;
}
}
int main(){try{
    const double e1=oscillator(100),e2=oscillator(200),e3=oscillator(400);
    require(e1/e2>3.8&&e1/e2<4.2&&e2/e3>3.8&&e2/e3<4.2,"Second-order convergence");
    // Explicitly time-dependent forcing, zero stiffness: a=t, u=t^3/6, v=t^2/2.
    ExplicitState timeState{};timeState.displacement={0};timeState.velocity={0};
    const InternalForceFunction noSpring=[](const DynamicVector&){return DynamicVector{0};};
    const ExternalForceFunction ramp=[](double t){return DynamicVector{t};};
    InitializeExplicitAcceleration(timeState,{1},{false},noSpring,ramp);
    for(int i=0;i<100;++i)advanceCentralDifference(timeState,{1},{false},.01,noSpring,ramp);
    near(timeState.acceleration[0],1,1e-12,"Force evaluated at new time");
    near(timeState.velocity[0],.5,1e-12,"Ramp velocity");
    near(timeState.displacement[0],1.0/6,2e-5,"Ramp displacement");
    // Fixed DOF with zero mass must not be divided by zero.
    ExplicitState fixed{};fixed.displacement={0,.01};fixed.velocity={0,0};
    const InternalForceFunction two=[](const DynamicVector& u){return DynamicVector{123,100*u[1]};};
    const ExternalForceFunction twoZero=[](double){return DynamicVector{0,0};};
    InitializeExplicitAcceleration(fixed,{0,1},{true,false},two,twoZero);
    for(int i=0;i<20;++i)advanceCentralDifference(fixed,{0,1},{true,false},.01,two,twoZero);
    require(fixed.displacement[0]==0&&fixed.velocity[0]==0&&fixed.acceleration[0]==0,"Fixed DOF drift");
    auto s=initial();InitializeExplicitAcceleration(s,{1},{false},spring,zero);const auto saved=s;
    for(double dt:{0.0,-.01,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
        rejects([&]{advanceCentralDifference(s,{1},{false},dt,spring,zero);});unchanged(s,saved);
    }
    rejects([&]{advanceCentralDifference(s,{0},{false},.01,spring,zero);});unchanged(s,saved);
    rejects([&]{advanceCentralDifference(s,{1},{},.01,spring,zero);});unchanged(s,saved);
    const ExternalForceFunction wrong=[](double){return DynamicVector{};};
    rejects([&]{advanceCentralDifference(s,{1},{false},.01,spring,wrong);});unchanged(s,saved);
    const ExternalForceFunction bad=[](double){return DynamicVector{std::numeric_limits<double>::infinity()};};
    rejects([&]{advanceCentralDifference(s,{1},{false},.01,spring,bad);});unchanged(s,saved);
    rejects([&]{InitializeExplicitAcceleration(s,{1},{false},spring,wrong);});unchanged(s,saved);
    std::cout<<"PASS: analytic oscillator, second-order convergence, energy, time-dependent force, constraints and invalid-input state protection.\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}return 0;}