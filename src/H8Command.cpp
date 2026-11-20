#include "nonlinearfem/H8Command.h"
#include "nonlinearfem/ModelIO_3D.h"
#include <iostream>
#include <stdexcept>
#include <string>
namespace NonlinearFEM {
int runH8Command(int argc,char** argv) {
    if(argc<2)return -1;
    const std::string command=argv[1];
    if(command!="--h8"&&command!="--dynamic-h8"&&command!="--help-h8")return -1;
    const bool dynamic=command=="--dynamic-h8";
    if(command=="--help-h8"||argc!=(dynamic?5:4)){
        std::cout<<"H8 static: NonlinearFEM --h8 model.nfem output_directory\n"
                   "H8 dynamic: NonlinearFEM --dynamic-h8 model.nfem dynamic.cfg output_directory\n";
        return command=="--help-h8"?0:1;
    }
    try {
        const std::filesystem::path directory=argv[dynamic?4:3];
        if(std::filesystem::exists(directory))throw std::runtime_error("Output directory already exists.");
        const auto model=readModelFile_3D(argv[2]);
        if(dynamic)runDynamicModel_3D(model,argv[3],directory);
        else {
            const auto result=solveLoadStepping_3D(model.mesh,model.material,model.loads,model.constraints,
                GlobalDisplacement_3D(3*model.mesh.nodes.size(),0),model.options);
            if(!result.converged)throw std::runtime_error("H8 analysis failed: "+result.failureReason);
            writeResults_3D(model,result,directory);
        }
        std::filesystem::copy_file(argv[2],directory/"model.nfem");
        std::cout<<"H8 output: "<<directory.string()<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return 1;}
}
}
