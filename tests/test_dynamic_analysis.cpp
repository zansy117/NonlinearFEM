#include "nonlinearfem/DynamicAnalysis.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace NonlinearFEM;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void near(double a,double b){require(std::isfinite(a)&&std::abs(a-b)<1e-12*std::max(1.0,std::abs(b)),"Numerical mismatch");}
void equalState(const ExplicitState& a,const ExplicitState& b){near(a.time,b.time);require(a.displacement.size()==b.displacement.size()&&a.velocity.size()==b.velocity.size()&&a.acceleration.size()==b.acceleration.size(),"State dimensions");for(std::size_t i=0;i<a.displacement.size();++i){near(a.displacement[i],b.displacement[i]);near(a.velocity[i],b.velocity[i]);near(a.acceleration[i],b.acceleration[i]);}}
template<class F>void reject(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,"Invalid input accepted");}
}
int main(){try{
 const Q4Mesh mesh={{{0,0},{1,0},{1,1},{0,1}},{{0,1,2,3}}};
 MaterialParameters mat{};mat.youngModulus=1000;mat.poissonRatio=.3;
 const std::vector<bool> fixed={true,true,false,false,false,false,true,true};
 ExplicitState initial{};initial.time=.2;initial.displacement.assign(8,0);initial.velocity.assign(8,0);initial.displacement[2]=initial.displacement[4]=1e-4;
 initial.acceleration.assign(8,999);const auto saved=initial;
 const auto mass=assembleGlobalLumpedMass(mesh,1,.2);
 const InternalForceFunction internal=[&](const DynamicVector& u){return assembleInternalForce(mesh,u,mat,.2);};
 const ExternalForceFunction zero=[](double){return DynamicVector(8,0);};
 const ExternalForceFunction ramp=[](double t){DynamicVector f(8,0);f[2]=f[4]=.01*t;return f;};
 for(const auto& external:{zero,ramp})for(std::size_t interval:{1u,10u,1000u}){
   DynamicAnalysisOptions options{};options.timeStep=1e-4;options.numberOfSteps=103;options.outputEvery=interval;
   const auto result=solveDynamicAnalysis(mesh,mat,1,.2,fixed,initial,external,options);
   equalState(initial,saved);
   ExplicitState direct=initial;InitializeExplicitAcceleration(direct,mass,fixed,internal,external);
   const std::size_t count=1+options.numberOfSteps/interval+(options.numberOfSteps%interval!=0);
   require(result.history.size()==count,"History count including final record");
   std::size_t recordIndex=0;
   const double initialEnergy=assembleGlobalEnergy(mesh,initial.displacement,mat,.2);
   double energyError=0;
   for(std::size_t step=0;step<=options.numberOfSteps;++step){
     if(step==0||step%interval==0||step==options.numberOfSteps){
       const auto& record=result.history.at(recordIndex++);equalState(record.state,direct);
       double kinetic=0;for(std::size_t i=0;i<8;++i){kinetic+=.5*mass[i]*direct.velocity[i]*direct.velocity[i];if(fixed[i])require(record.state.displacement[i]==0&&record.state.velocity[i]==0&&record.state.acceleration[i]==0,"Constraint drift");}
       near(record.kineticEnergy,kinetic);near(record.strainEnergy,assembleGlobalEnergy(mesh,direct.displacement,mat,.2));
       energyError=std::max(energyError,std::abs(record.kineticEnergy+record.strainEnergy-initialEnergy)/initialEnergy);
     }
     if(step<options.numberOfSteps)advanceCentralDifference(direct,mass,fixed,options.timeStep,internal,external);
   }
   equalState(result.finalState,direct);equalState(result.history.back().state,result.finalState);
   near(result.finalState.time,.2103);
   if(external.target_type()==zero.target_type())require(energyError<1e-4,"Free-vibration energy error");
 }
 DynamicAnalysisOptions options{};
 for(double dt:{0.0,-1.0,std::numeric_limits<double>::infinity()}){options.timeStep=dt;reject([&]{(void)solveDynamicAnalysis(mesh,mat,1,.2,fixed,initial,zero,options);});}
 options={};options.numberOfSteps=0;reject([&]{(void)solveDynamicAnalysis(mesh,mat,1,.2,fixed,initial,zero,options);});
 options={};options.outputEvery=0;reject([&]{(void)solveDynamicAnalysis(mesh,mat,1,.2,fixed,initial,zero,options);});
 options={};reject([&]{(void)solveDynamicAnalysis(mesh,mat,0,.2,fixed,initial,zero,options);});
 reject([&]{(void)solveDynamicAnalysis(mesh,mat,1,.2,{},initial,zero,options);});
 reject([&]{(void)solveDynamicAnalysis(mesh,mat,1,.2,fixed,initial,ExternalForceFunction{},options);});
 equalState(initial,saved);
 std::cout<<"PASS: direct-step equivalence, nonzero initial time, ramp load, initialization, record scheduling, energies, constraints, input immutability and invalid controls.\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}return 0;}