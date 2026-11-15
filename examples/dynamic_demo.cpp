#include "nonlinearfem/DynamicOutput.h"
#include <iostream>
#include <stdexcept>
using namespace NonlinearFEM;
int main(int argc,char** argv) {
    if(argc!=2){std::cerr<<"Usage: dynamic_demo new_output_directory\n";return 1;}
    try {
        if(std::filesystem::exists(argv[1]))throw std::runtime_error("Output directory already exists.");
        Q4Mesh mesh;
        for(int j=0;j<=2;++j)for(int i=0;i<=4;++i)mesh.nodes.push_back({double(i),.5*j});
        for(std::size_t j=0;j<2;++j)for(std::size_t i=0;i<4;++i){const auto a=5*j+i;mesh.elements.push_back({a,a+1,a+6,a+5});}
        MaterialParameters material{};material.youngModulus=1000;material.poissonRatio=.3;
        ExplicitState initial{};initial.displacement.assign(30,0);initial.velocity.assign(30,0);
        std::vector<bool> fixed(30,false);
        for(std::size_t n=0;n<15;++n){
            const double x=mesh.nodes[n][0],y=mesh.nodes[n][1]-.5;
            if(n%5==0)fixed[2*n]=fixed[2*n+1]=true;
            initial.displacement[2*n]=-y*2e-4*x/16;
            initial.displacement[2*n+1]=1e-4*x*x/16;
        }
        const ExternalForceFunction external=[](double){return DynamicVector(30,0);};
        DynamicAnalysisOptions options{};options.timeStep=2e-4;options.numberOfSteps=10000;options.outputEvery=25;
        const auto result=solveDynamicAnalysis(mesh,material,1,.2,fixed,initial,external,options);
        writeDynamicResults(mesh,result,9,argv[1]);
        std::cout<<"Exported "<<result.history.size()<<" frames; observer=node 9; dt=0.0002; duration=2.\n";
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}
