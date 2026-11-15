#include "nonlinearfem/Assembly.h"
#include "nonlinearfem/ExplicitDynamics.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace NonlinearFEM;
namespace {
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
void near(double a,double b){require(std::isfinite(a)&&std::abs(a-b)<1e-12*std::max(1.0,std::abs(b)),"Analytical time-step mismatch");}
template<class F>void reject(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,"Invalid input accepted");}
Q4Mesh mesh(int nx,int ny){Q4Mesh m;for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i)m.nodes.push_back({4.0*i/nx,1.0*j/ny});for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){auto a=static_cast<std::size_t>(j*(nx+1)+i);auto b=a+static_cast<std::size_t>(nx+1);m.elements.push_back({a,a+1,b+1,b});}return m;}
MaterialParameters material(){MaterialParameters m{};m.youngModulus=1000;m.poissonRatio=.3;return m;}
std::vector<bool> constraints(const Q4Mesh& m){std::vector<bool> c(2*m.nodes.size(),false);for(std::size_t i=0;i<m.nodes.size();++i)if(m.nodes[i][0]==0)c[2*i]=c[2*i+1]=true;return c;}
double estimate(const Q4Mesh& m){const auto mat=material();return estimateStableTimeStep(assembleGlobalSystem(m,GlobalDisplacement(2*m.nodes.size(),0),mat,.2).tangentStiffness,assembleGlobalLumpedMass(m,1,.2),constraints(m),.8);}
double vibration(const Q4Mesh& m,double dt,int steps){
 const auto mat=material();const auto mass=assembleGlobalLumpedMass(m,1,.2);const auto fixed=constraints(m);const auto N=mass.size();
 ExplicitState s{};s.displacement.assign(N,0);s.velocity.assign(N,0);
 for(std::size_t n=0;n<m.nodes.size();++n){double x=m.nodes[n][0],y=m.nodes[n][1]-.5;s.displacement[2*n]=-y*2e-6*x/16;s.displacement[2*n+1]=1e-6*x*x/16;}
 const InternalForceFunction internal=[&](const DynamicVector& u){return assembleInternalForce(m,u,mat,.2);};const ExternalForceFunction external=[N](double){return DynamicVector(N,0);};
 const auto energy=[&](){double t=0;for(std::size_t i=0;i<N;++i)t+=.5*mass[i]*s.velocity[i]*s.velocity[i];return t+assembleGlobalEnergy(m,s.displacement,mat,.2);};
 InitializeExplicitAcceleration(s,mass,fixed,internal,external);const double initial=energy();require(initial>0,"Zero initial energy");double error=0;
 for(int step=0;step<steps;++step){advanceCentralDifference(s,mass,fixed,dt,internal,external);for(std::size_t i=0;i<N;++i){require(std::isfinite(s.displacement[i])&&std::abs(s.displacement[i])<1e-5,"Unbounded Q4 vibration");if(fixed[i])require(s.displacement[i]==0&&s.velocity[i]==0&&s.acceleration[i]==0,"Fixed DOF drift");}const double e=energy();require(std::isfinite(e),"Non-finite energy");error=std::max(error,std::abs(e-initial)/initial);}
 require(error<.1,"Excessive Q4 energy error");return error;
}
}
int main(){try{
 near(estimateStableTimeStep({{100}},{1},{false},.8),.16);
 // Nonuniform masses and coupled stiffness: normalized matrix [[4,-2],[-2,4]], eigenvalues 2,6.
 near(estimateStableTimeStep({{4,-4},{-4,16}},{1,4},{false,false},.8),1.6/std::sqrt(6.0));
 near(estimateStableTimeStep({{100}},{4},{false},.8),.32);
 near(estimateStableTimeStep({{400}},{1},{false},.8),.08);
 near(estimateStableTimeStep({{100}},{1},{false},.4),.08);
 near(estimateStableTimeStep({{1e12,50},{50,100}},{0,1},{true,false},.8),.16);
 const double nan=std::numeric_limits<double>::quiet_NaN();
 reject([]{(void)estimateStableTimeStep({}, {}, {}, .8);});
 reject([]{(void)estimateStableTimeStep({{1}}, {}, {false}, .8);});
 reject([]{(void)estimateStableTimeStep({{}},{1},{false},.8);});
 reject([]{(void)estimateStableTimeStep({{1}},{0},{false},.8);});
 reject([]{(void)estimateStableTimeStep({{1}},{1},{true},.8);});
 reject([]{(void)estimateStableTimeStep({{0}},{1},{false},.8);});
 reject([]{(void)estimateStableTimeStep({{1,2},{0,1}},{1,1},{false,false},.8);});
 for(double factor:{0.0,1.0,-1.0,nan})reject([&]{(void)estimateStableTimeStep({{1}},{1},{false},factor);});
 reject([&]{(void)estimateStableTimeStep({{nan}},{1},{false},.8);});
 reject([&]{(void)estimateStableTimeStep({{1}},{nan},{false},.8);});
 const auto coarse=mesh(4,2),fine=mesh(8,4);const double dt=estimate(coarse),fineDt=estimate(fine);
 std::cout<<"Q4 estimated steps: 4x2="<<dt<<", 8x4="<<fineDt<<'\n';require(fineDt>0&&fineDt<dt,"Refined mesh did not reduce estimated step");
 const double a=vibration(coarse,dt,200),b=vibration(coarse,dt/2,400),c=vibration(coarse,dt/4,800);
 std::cout<<"Q4 maximum relative energy errors: "<<a<<", "<<b<<", "<<c<<'\n';
 require(a/b>3&&a/b<5&&b/c>3&&b/c<5,"Estimated-step energy refinement failed");
 std::cout<<"PASS: analytical bounds, mass/stiffness scaling, constraints, invalid inputs, mesh refinement and Q4 dynamics.\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}return 0;}