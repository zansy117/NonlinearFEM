#include "nonlinearfem/ModelIO_3D.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace NonlinearFEM {
void runDynamicModel_3D(const AnalysisModel_3D& model,
    const std::filesystem::path& configuration,
    const std::filesystem::path& output) {
    std::ifstream input(configuration);
    if(!input)throw std::invalid_argument("Cannot open dynamic configuration.");
    std::string magic;int version;
    if(!(input>>magic>>version)||magic!="NFEM_DYNAMIC"||version!=1)
        throw std::invalid_argument("Expected NFEM_DYNAMIC 1 header.");
    std::map<std::string,std::string> values;
    std::string key,value;
    while(input>>key){
        if(!(input>>value)||!values.emplace(key,value).second)
            throw std::invalid_argument("Invalid or duplicate dynamic option.");
    }
    const std::vector<std::string> keys={"density","time_step","duration","output_every",
        "initial_amplitude","observation_node","load","initial_mode"};
    if(values.size()!=keys.size())throw std::invalid_argument("Incorrect dynamic option count.");
    for(const auto& name:keys)if(!values.count(name))throw std::invalid_argument("Missing dynamic option: "+name);
    const auto number=[&](const std::string& name){
        std::size_t used=0;double n=std::stod(values.at(name),&used);
        if(used!=values.at(name).size()||!std::isfinite(n))throw std::invalid_argument("Invalid number: "+name);
        return n;
    };
    const double density=number("density"),requestedDt=number("time_step"),duration=number("duration"),amplitude=number("initial_amplitude");
    const double every=number("output_every"),observer=number("observation_node");
    if(density<=0||requestedDt<0||duration<=0||every<1||every!=std::floor(every)||every>1e6||
        observer<0||observer!=std::floor(observer)||observer>=model.mesh.nodes.size())
        throw std::invalid_argument("Invalid dynamic density, step, duration, output interval or observer.");
    const auto load=values.at("load"),mode=values.at("initial_mode");
    if((load!="free"&&load!="ramp"&&load!="step")||(mode!="bending"&&mode!="axial"&&mode!="shear"&&mode!="sine_bending"))
        throw std::invalid_argument("Unknown loading or perturbation mode.");
    const auto N=3*model.mesh.nodes.size();
    std::vector<bool> fixed(N,false);
    for(const auto& c:model.constraints){
        if(c.value!=0)throw std::invalid_argument("Dynamic analysis supports only fixed zero constraints.");
        fixed[3*c.node+c.direction]=true;
    }
    ExplicitState initial{};initial.displacement.assign(N,0);initial.velocity.assign(N,0);
    double xmin=model.mesh.nodes[0][0],xmax=xmin,ymin=model.mesh.nodes[0][1],ymax=ymin;
    for(const auto& p:model.mesh.nodes){xmin=std::min(xmin,p[0]);xmax=std::max(xmax,p[0]);ymin=std::min(ymin,p[1]);ymax=std::max(ymax,p[1]);}
    const double length=xmax-xmin;
    if(length<=0)throw std::invalid_argument("Initial perturbation requires nonzero X extent.");
    for(std::size_t n=0;n<model.mesh.nodes.size();++n){
        const double x=model.mesh.nodes[n][0]-xmin,y=model.mesh.nodes[n][1]-(ymin+ymax)/2;
        initial.displacement[3*n]=mode=="axial"?amplitude*x/length:-y*2*amplitude*x/(length*length);
        initial.displacement[3*n+1]=mode=="axial"?0:amplitude*x*x/(length*length);
        if(mode=="shear"){
            const double height=ymax-ymin;
            if(height<=0)throw std::invalid_argument("Shear perturbation requires nonzero Y extent.");
            initial.displacement[3*n]=amplitude*(model.mesh.nodes[n][1]-ymin)/height;
            initial.displacement[3*n+1]=0;
        } else if(mode=="sine_bending"){
            const double pi=std::acos(-1.0);
            initial.displacement[3*n]=-y*amplitude*pi/length*std::cos(pi*x/length);
            initial.displacement[3*n+1]=amplitude*std::sin(pi*x/length);
        }
        for(std::size_t d=0;d<3;++d)if(fixed[3*n+d])initial.displacement[3*n+d]=0;
    }
    const auto mass=assembleGlobalLumpedMass_3D(model.mesh,density);
    const auto system=assembleGlobalSystem_3D(model.mesh,initial.displacement,model.material);
    // Remove only roundoff-level asymmetry in the estimate's matrix copy.
    auto estimateStiffness=system.tangentStiffness;
    double stiffnessScale=1.0;
    for(const auto& row:estimateStiffness)for(double k:row)stiffnessScale=std::max(stiffnessScale,std::abs(k));
    for(std::size_t i=0;i<N;++i)for(std::size_t j=i+1;j<N;++j){
        const double a=estimateStiffness[i][j],b=estimateStiffness[j][i];
        if(std::abs(a/stiffnessScale-b/stiffnessScale)>1e-12)
            throw std::runtime_error("Initial stiffness has significant asymmetry.");
        estimateStiffness[i][j]=estimateStiffness[j][i]=.5*a+.5*b;
    }
    const double bound=estimateStableTimeStep(estimateStiffness,mass,fixed,.8);
    const double selected=requestedDt==0?bound*.2:requestedDt;
    if(selected>bound)throw std::invalid_argument("Requested time step exceeds initial-state conservative estimate; reduce time_step or use 0 for automatic selection.");
    const double count=std::ceil(duration/selected);
    if(!std::isfinite(count)||count<1||count>1000000)throw std::invalid_argument("Too many time steps (maximum 1000000). Reduce duration or model stiffness.");
    DynamicAnalysisOptions options{};options.numberOfSteps=static_cast<std::size_t>(count);
    options.timeStep=duration/options.numberOfSteps;options.outputEvery=static_cast<std::size_t>(every);
    const auto records=2+options.numberOfSteps/options.outputEvery;
    if(records>2000000/model.mesh.nodes.size())throw std::invalid_argument("Too many saved node states. Increase output_every.");
    const auto force=assembleExternalForce_3D(model.mesh.nodes.size(),model.loads);
    const ExternalForceFunction external=[force,load,duration](double t){
        auto result=force;
        const double factor=load=="free"?0:load=="step"?1:std::min(1.0,std::max(0.0,t/duration));
        for(auto& f:result)f*=factor;
        return result;
    };
    const auto result=solveDynamicAnalysis_3D(model.mesh,model.material,density,fixed,initial,external,options);
    writeDynamicResults_3D(model.mesh,result,static_cast<std::size_t>(observer),output);
    std::filesystem::copy_file(configuration,output/"dynamic.cfg");
    std::ofstream report(output/"analysis.txt");report.exceptions(std::ios::failbit|std::ios::badbit);
    report<<std::setprecision(17)<<"analysis dynamic\ndimension 3\nload "<<load<<"\ndensity "<<density<<"\ntime_step "<<options.timeStep
        <<"\ninitial_step_estimate "<<bound<<"\nnumber_of_steps "<<options.numberOfSteps<<"\nduration "<<duration
        <<"\nobservation_node "<<observer<<"\n";
    std::cout<<"Dynamic analysis completed: dt="<<options.timeStep<<", steps="<<options.numberOfSteps
        <<", saved frames="<<result.history.size()<<". Initial-state estimate is not a nonlinear stability guarantee.\n";
}
}
