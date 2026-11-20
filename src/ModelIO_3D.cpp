#include "nonlinearfem/ModelIO_3D.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace NonlinearFEM {
namespace {
struct Reader {
    std::istream& input;
    std::size_t line = 0;
    [[noreturn]] void fail(const std::string& message) const {
        throw std::runtime_error("Model line " + std::to_string(line) + ": " + message);
    }
    bool next(std::istringstream& row) {
        std::string text;
        while (std::getline(input, text)) {
            ++line;
            text = text.substr(0, text.find('#'));
            row.clear(); row.str(text);
            row >> std::ws;
            if (!row.eof()) return true;
        }
        if (input.bad()) fail("input read failed");
        return false;
    }
    void end(std::istringstream& row) const {
        row >> std::ws;
        if (!row.eof()) fail("unexpected extra value");
    }
    template<class... T> void values(std::istringstream& row, T&... values) const {
        if (!(row >> ... >> values)) fail("missing or invalid value");
        end(row);
    }
    void record(std::istringstream& row) {
        if (!next(row)) fail("unexpected end of file");
    }
};
void finite(double value, const Reader& r) {
    if (!std::isfinite(value)) r.fail("all numeric values must be finite");
}
std::ofstream output(const std::filesystem::path& path) {
    std::ofstream stream(path);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << std::setprecision(17);
    return stream;
}
}

AnalysisModel_3D readModel_3D(std::istream& input) {
    Reader r{input};
    AnalysisModel_3D m{};
    std::istringstream row;
    r.record(row);
    std::string magic; int version = 0;
    r.values(row, magic, version);
    if (magic != "NFEM3D" || version != 1) r.fail("expected NFEM3D 1 header");
    std::set<std::string> seen;
    bool ended = false;
    while (r.next(row)) {
        std::string key; row >> key;
        if (key == "end") { r.end(row); ended = true; break; }
        if (!seen.insert(key).second) r.fail("duplicate section: " + key);
        if (key == "material") {
            r.values(row, m.material.youngModulus, m.material.poissonRatio);
            finite(m.material.youngModulus, r); finite(m.material.poissonRatio, r);
            if (m.material.youngModulus <= 0 || m.material.poissonRatio <= -1 ||
                m.material.poissonRatio >= 0.5) r.fail("invalid material parameters");
        } else if (key == "newton") {
            long long iterations;
            r.values(row, iterations, m.options.newtonOptions.residualTolerance,
                     m.options.newtonOptions.displacementTolerance);
            finite(m.options.newtonOptions.residualTolerance, r);
            finite(m.options.newtonOptions.displacementTolerance, r);
            if (iterations < 1 || iterations > 10000 ||
                m.options.newtonOptions.residualTolerance <= 0 ||
                m.options.newtonOptions.displacementTolerance <= 0) r.fail("invalid Newton options");
            m.options.newtonOptions.maxIterations = static_cast<std::size_t>(iterations);
        } else if (key == "stepping") {
            auto& o = m.options;
            long long retries, steps;
            r.values(row, o.initialIncrement, o.minimumIncrement, o.maximumIncrement,
                     o.cutbackFactor, o.growthFactor, retries, steps);
            for (double v : {o.initialIncrement,o.minimumIncrement,o.maximumIncrement,
                             o.cutbackFactor,o.growthFactor}) finite(v,r);
            if (!(o.minimumIncrement > 0 && o.initialIncrement >= o.minimumIncrement &&
                  o.maximumIncrement >= o.initialIncrement && o.maximumIncrement <= 1 &&
                  o.cutbackFactor > 0 && o.cutbackFactor < 1 && o.growthFactor >= 1 &&
                  retries >= 0 && retries <= 10000 && steps > 0 && steps <= 100000))
                r.fail("invalid stepping options");
            o.maxRetries=static_cast<std::size_t>(retries);
            o.maxLoadSteps=static_cast<std::size_t>(steps);
        } else if (key == "nodes" || key == "elements" || key == "loads" || key == "constraints") {
            long long count;
            r.values(row,count);
            if (count < 0 || count > 1000000) r.fail("invalid record count");
            for (long long i=0; i<count; ++i) {
                r.record(row);
                if (key == "nodes") {
                    H8NodeCoordinate point{};
                    r.values(row,point[0],point[1],point[2]);
                    finite(point[0],r); finite(point[1],r); finite(point[2],r);
                    m.mesh.nodes.push_back(point);
                } else if (key == "elements") {
                    std::array<long long,8> indices{};
                    r.values(row,indices[0],indices[1],indices[2],indices[3],indices[4],indices[5],indices[6],indices[7]);
                    H8ElementConnectivity element{};
                    for (std::size_t a=0; a<8; ++a) {
                        if (indices[a]<0) r.fail("negative node index");
                        element[a]=static_cast<std::size_t>(indices[a]);
                    }
                    m.mesh.elements.push_back(element);
                } else {
                    long long node, direction; double value;
                    r.values(row,node,direction,value); finite(value,r);
                    if (node<0 || direction<0 || direction>2) r.fail("invalid node or direction (use 0=x, 1=y, 2=z)");
                    const auto n=static_cast<std::size_t>(node);
                    const auto d=static_cast<std::size_t>(direction);
                    if (key == "loads") m.loads.push_back({n,d,value});
                    else m.constraints.push_back({n,d,value});
                }
            }
        } else r.fail("unknown section: " + key);
    }
    if (!ended) r.fail("missing end marker");
    if (r.next(row)) r.fail("content after end marker");
    for (const auto* key : {"nodes","elements","material","loads","constraints"})
        if (!seen.count(key)) r.fail(std::string("missing section: ")+key);
    // Dense solver: protect the CLI against inadvertently huge matrix allocations.
    if (m.mesh.nodes.size()>600) r.fail("dense solver supports at most 600 input nodes");
    std::set<std::size_t> constrained;
    for (const auto& c : m.constraints) {
        if (c.node>=m.mesh.nodes.size()) r.fail("constraint node out of range");
        if (!constrained.insert(3*c.node+c.direction).second) r.fail("duplicate constrained DOF");
    }
    for (const auto& load : m.loads)
        if (load.node>=m.mesh.nodes.size()) r.fail("load node out of range");
    // Reuse mesh/material/Jacobian validation before exposing the model to Solver.
    const auto initial = assembleInternalForce_3D(m.mesh,
        GlobalDisplacement_3D(3*m.mesh.nodes.size(),0.0),m.material);
    (void)initial;
    return m;
}
AnalysisModel_3D readModelFile_3D(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open model: "+path.string());
    return readModel_3D(input);
}


}
