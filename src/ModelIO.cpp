#include "nonlinearfem/ModelIO.h"
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

AnalysisModel readModel(std::istream& input) {
    Reader r{input};
    AnalysisModel m{};
    std::istringstream row;
    r.record(row);
    std::string magic; int version = 0;
    r.values(row, magic, version);
    if (magic != "NFEM" || version != 1) r.fail("expected NFEM 1 header");
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
        } else if (key == "thickness") {
            r.values(row, m.thickness); finite(m.thickness, r);
            if (m.thickness <= 0) r.fail("thickness must be positive");
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
                    Q4NodeCoordinate point{};
                    r.values(row,point[0],point[1]);
                    finite(point[0],r); finite(point[1],r);
                    m.mesh.nodes.push_back(point);
                } else if (key == "elements") {
                    std::array<long long,4> indices{};
                    r.values(row,indices[0],indices[1],indices[2],indices[3]);
                    Q4ElementConnectivity element{};
                    for (std::size_t a=0; a<4; ++a) {
                        if (indices[a]<0) r.fail("negative node index");
                        element[a]=static_cast<std::size_t>(indices[a]);
                    }
                    m.mesh.elements.push_back(element);
                } else {
                    long long node, direction; double value;
                    r.values(row,node,direction,value); finite(value,r);
                    if (node<0 || direction<0 || direction>1) r.fail("invalid node or direction (use 0=x, 1=y)");
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
    for (const auto* key : {"nodes","elements","material","thickness","loads","constraints"})
        if (!seen.count(key)) r.fail(std::string("missing section: ")+key);
    // Dense solver: protect the CLI against inadvertently huge matrix allocations.
    if (m.mesh.nodes.size()>2000) r.fail("dense solver supports at most 2000 input nodes");
    std::set<std::size_t> constrained;
    for (const auto& c : m.constraints) {
        if (c.node>=m.mesh.nodes.size()) r.fail("constraint node out of range");
        if (!constrained.insert(2*c.node+c.direction).second) r.fail("duplicate constrained DOF");
    }
    for (const auto& load : m.loads)
        if (load.node>=m.mesh.nodes.size()) r.fail("load node out of range");
    // Reuse mesh/material/Jacobian validation before exposing the model to Solver.
    const auto initial = assembleGlobalSystem(m.mesh,
        GlobalDisplacement(2*m.mesh.nodes.size(),0.0),m.material,m.thickness);
    (void)initial;
    return m;
}
AnalysisModel readModelFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open model: "+path.string());
    return readModel(input);
}

void writeResults(const AnalysisModel& m, const LoadSteppingResult& result,
                  const std::filesystem::path& directory) {
    if (!result.converged || std::abs(result.finalLoadFactor-1.0)>1e-12 ||
        result.displacement.size()!=2*m.mesh.nodes.size())
        throw std::runtime_error("Only fully converged results can be exported");
    for (double u : result.displacement)
        if (!std::isfinite(u)) throw std::runtime_error("Non-finite result");
    const auto assembly=assembleGlobalSystem(m.mesh,result.displacement,m.material,m.thickness);
    const auto external=assembleExternalForce(m.mesh.nodes.size(),m.loads);
    std::vector<bool> fixed(result.displacement.size(),false);
    for (const auto& c : m.constraints) fixed[2*c.node+c.direction]=true;
    if (!std::filesystem::create_directories(directory))
        throw std::runtime_error("Output directory already exists");
    {
        auto csv=output(directory/"nodes.csv");
        csv<<"node,x,y,ux,uy,displacement_magnitude,deformed_x,deformed_y,reaction_x,reaction_y\n";
        for (std::size_t n=0;n<m.mesh.nodes.size();++n) {
            const double ux=result.displacement[2*n], uy=result.displacement[2*n+1];
            csv<<n<<','<<m.mesh.nodes[n][0]<<','<<m.mesh.nodes[n][1]<<','<<ux<<','<<uy
               <<','<<std::hypot(ux,uy)<<','<<m.mesh.nodes[n][0]+ux<<','<<m.mesh.nodes[n][1]+uy
               <<','<<(fixed[2*n] ? assembly.internalForce[2*n]-external[2*n] : 0.0)
               <<','<<(fixed[2*n+1] ? assembly.internalForce[2*n+1]-external[2*n+1] : 0.0)<<'\n';
        }
    }
    {
        auto csv=output(directory/"history.csv");
        csv<<"step,load_factor,increment,newton_iterations,retries,relative_residual\n";
        for (std::size_t i=0;i<result.history.size();++i) {
            const auto& h=result.history[i];
            csv<<i+1<<','<<h.loadFactor<<','<<h.loadIncrement<<','<<h.iterations
               <<','<<h.retries<<','<<h.residualNorm<<'\n';
        }
    }
    auto html=output(directory/"results.html");
    html<<R"HTML(<!doctype html><html lang="zh-CN"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>NonlinearFEM · 计算结果</title>
<style>body{font:16px system-ui;margin:0;background:#f3f6fa;color:#172c40}main{max-width:1100px;margin:32px auto;padding:0 20px}h1{font-size:28px}.panel{background:white;padding:20px;border-radius:12px;margin:16px 0}label{display:inline-block;margin:8px 24px 8px 0}canvas{width:100%;height:540px;display:block}input[type=number]{width:90px}#legend{height:12px;background:linear-gradient(90deg,hsl(240,80%,50%),hsl(180,80%,50%),hsl(120,80%,50%),hsl(60,80%,50%),hsl(0,80%,50%))}.ends{display:flex;justify-content:space-between}a{color:#175ea8}#detail{font-family:monospace;min-height:42px}small{color:#526579}</style>
<main><h1>NonlinearFEM / 位移与变形</h1><p>最终载荷系数 1 · 已收敛 · 单位沿用输入模型</p>
<div class="panel"><label>位移分量 <select id="field"><option value="mag">位移大小 |u|</option><option value="ux">水平位移 ux</option><option value="uy">竖向位移 uy</option></select></label>
<label>变形放大倍数 <input id="scale" type="number" min="0" max="10000" value="1" step="0.5"></label>
<label><input id="original" type="checkbox" checked> 显示原始网格</label>
<canvas id="view" aria-label="原始网格与位移着色的变形网格"></canvas>
<div id="legend"></div><div class="ends"><span id="min"></span><span id="max"></span></div>
<p><small>灰色虚线：原始网格；彩色：变形网格。单元按节点位移平均值着色，图例范围采用同一平均值；不是应力云图。移动鼠标查看最近节点。</small></p><div id="detail"></div></div>
<div class="panel"><a href="nodes.csv" download>下载节点位移与支座反力</a> · <a href="history.csv" download>下载加载记录</a><p>仅导出最终状态。放大倍数只改变显示，不改变计算结果。</p></div></main><script>
const nodes=)HTML";
    html<<'[';
    for (std::size_t n=0;n<m.mesh.nodes.size();++n) {
        if(n) html<<',';
        html<<'['<<m.mesh.nodes[n][0]<<','<<m.mesh.nodes[n][1]<<','
            <<result.displacement[2*n]<<','<<result.displacement[2*n+1]<<']';
    }
    html<<"], elements=[";
    for (std::size_t i=0;i<m.mesh.elements.size();++i) {
        if(i) html<<',';
        const auto& e=m.mesh.elements[i];
        html<<'['<<e[0]<<','<<e[1]<<','<<e[2]<<','<<e[3]<<']';
    }
    html<<R"HTML(];
const canvas=document.getElementById('view'),ctx=canvas.getContext('2d');
const field=document.getElementById('field'),scale=document.getElementById('scale'),original=document.getElementById('original');
let screen=[];
function draw(){
 const w=canvas.clientWidth,h=540,dpr=window.devicePixelRatio||1;
 canvas.width=w*dpr;canvas.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
 const factor=Number(scale.value);if(!Number.isFinite(factor)||factor<0||factor>10000)return;
 const points=nodes.map(n=>[n[0]+factor*n[2],n[1]+factor*n[3]]);
 let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;
 for(let i=0;i<nodes.length;i++)for(const p of (original.checked?[points[i],nodes[i]]:[points[i]])){
  xmin=Math.min(xmin,p[0]);xmax=Math.max(xmax,p[0]);ymin=Math.min(ymin,p[1]);ymax=Math.max(ymax,p[1]);}
 const zoom=Math.min((w-60)/Math.max(xmax-xmin,1e-12),(h-60)/Math.max(ymax-ymin,1e-12));
 const project=p=>[(w-(xmax-xmin)*zoom)/2+(p[0]-xmin)*zoom,(h+(ymax-ymin)*zoom)/2-(p[1]-ymin)*zoom];
 const value=n=>field.value==='ux'?n[2]:field.value==='uy'?n[3]:Math.hypot(n[2],n[3]);
 const values=elements.map(e=>e.reduce((sum,i)=>sum+value(nodes[i]),0)/4);
 const lo=Math.min(...values),hi=Math.max(...values);
 document.getElementById('min').textContent=lo.toExponential(5);document.getElementById('max').textContent=hi.toExponential(5);
 const polygon=(e,ps)=>{ctx.beginPath();e.forEach((n,k)=>{const p=project(ps[n]);k?ctx.lineTo(...p):ctx.moveTo(...p)});ctx.closePath()};
 elements.forEach((e,k)=>{polygon(e,points);const t=hi===lo?0.5:(values[k]-lo)/(hi-lo);ctx.fillStyle=`hsl(${240*(1-t)},80%,50%)`;ctx.fill();ctx.strokeStyle='#17314b';ctx.lineWidth=0.7;ctx.stroke()});
 if(original.checked){ctx.setLineDash([4,4]);ctx.strokeStyle='#6a7785';ctx.lineWidth=1;elements.forEach(e=>{polygon(e,nodes);ctx.stroke()});ctx.setLineDash([])}
 screen=points.map(project);
}
canvas.addEventListener('mousemove',event=>{const rect=canvas.getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top;let best=-1,distance=Infinity;screen.forEach((p,i)=>{const d=Math.hypot(p[0]-x,p[1]-y);if(d<distance){distance=d;best=i}});document.getElementById('detail').textContent=best<0?'':`节点 ${best} | X=(${nodes[best][0]}, ${nodes[best][1]}) | ux=${nodes[best][2].toExponential(6)} | uy=${nodes[best][3].toExponential(6)}`;});
field.addEventListener('change',draw);scale.addEventListener('input',draw);original.addEventListener('change',draw);window.addEventListener('resize',draw);draw();
</script></html>)HTML";
}
}
