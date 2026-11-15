#include "nonlinearfem/ModelIO.h"
#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace NonlinearFEM;
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void reject(const std::string& text) {
    std::istringstream input(text);
    bool rejected=false;
    try { (void)readModel(input); } catch (const std::exception&) { rejected=true; }
    require(rejected,"Invalid model accepted");
}
std::string contents(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input),"Missing exported file");
    return {std::istreambuf_iterator<char>(input),std::istreambuf_iterator<char>()};
}
}
int main(int argc,char** argv) {
    try {
        require(argc==2,"Expected model fixture path");
        const auto fixture=contents(argv[1]);
        auto model=readModelFile(argv[1]);
        auto result=solveLoadStepping(model.mesh,model.material,model.loads,model.constraints,
            GlobalDisplacement(18,0.0),model.options,model.thickness);
        require(result.converged,"File-based tension solve failed");
        const double transverse=std::sqrt(1-0.3*(1.2*1.2-1));
        for (std::size_t n=0;n<9;++n) {
            require(std::abs(result.displacement[2*n]-0.2*model.mesh.nodes[n][0])<1e-8,
                    "Imported model axial solution mismatch");
            require(std::abs(result.displacement[2*n+1]-(transverse-1)*model.mesh.nodes[n][1])<1e-8,
                    "Imported model transverse solution mismatch");
        }
        for (const auto& replacement : {std::string("99 0 0"),std::string("0 2 0"),std::string("-1 0 0")}) {
            auto invalid=fixture;
            invalid.replace(invalid.find("0 0 0\n"),6,replacement);
            reject(invalid);
        }
        reject(fixture+"garbage\n");
        reject(fixture.substr(0,fixture.find("end\n")));
        auto duplicate=fixture;
        duplicate.replace(duplicate.find("6 0 0"),5,"0 0 0");
        reject(duplicate);
        auto badElement=fixture;
        badElement.replace(badElement.find("0 1 4 3"),7,"0 1 4 99");
        reject(badElement);
        const auto directory=std::filesystem::temp_directory_path()/
            ("nfem_io_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
        writeResults(model,result,directory);
        const auto csv=contents(directory/"nodes.csv");
        std::istringstream rows(csv); std::string line; std::getline(rows,line);
        double rx=0,ry=0; std::size_t count=0;
        while (std::getline(rows,line)) {
            for(char& c:line) if(c==',') c=' ';
            std::istringstream row(line);
            std::size_t id; double x,y,ux,uy,mag,dx,dy,fx,fy;
            require(static_cast<bool>(row>>id>>x>>y>>ux>>uy>>mag>>dx>>dy>>fx>>fy),"Invalid CSV row");
            require(id==count && std::abs(ux-result.displacement[2*id])<1e-14,"CSV displacement mismatch");
            require(std::abs(dx-x-ux)<1e-14 && std::abs(dy-y-uy)<1e-14,"CSV deformed position mismatch");
            rx+=fx;ry+=fy;++count;
        }
        require(count==9 && std::abs(rx+11088)<1e-4 && std::abs(ry)<1e-4,"CSV reaction balance failed");
        require(contents(directory/"history.csv").find("load_factor")!=std::string::npos,"Missing history");
        require(contents(directory/"results.html").find("const nodes=[")!=std::string::npos,"Missing viewer data");
        bool protectedOutput=false;
        try { writeResults(model,result,directory); } catch(const std::exception&) { protectedOutput=true; }
        require(protectedOutput,"Existing output overwritten");
        // Remove only the three files created by this test, never a recursive directory deletion.
        for(const auto* name:{"nodes.csv","history.csv","results.html"}) std::filesystem::remove(directory/name);
        std::filesystem::remove(directory);
        std::cout<<"Model parsing, analytic solve, CSV reactions and output protection passed.\n";
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
    return 0;
}
