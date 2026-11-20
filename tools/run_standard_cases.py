"""Run the manual's frozen reference cases and save auditable results.

Run from any directory: python tools/run_standard_cases.py --build-dir build
This reuses existing solver/tests and does not alter them or overwrite a report.
"""
import argparse,csv,hashlib,json,math,platform,subprocess,sys
from datetime import datetime,timezone,timedelta
from pathlib import Path

def main(default_project=None):
    parser=argparse.ArgumentParser()
    parser.add_argument('--project',type=Path,default=default_project or Path(__file__).resolve().parents[1])
    parser.add_argument('--build-dir',type=Path,default=Path('build'))
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();root=args.project.resolve();build=args.build_dir if args.build_dir.is_absolute() else root/args.build_dir
    now=datetime.now(timezone(timedelta(hours=8)))
    output=(args.output or root/'results'/('standard_'+now.strftime('%Y%m%d_%H%M%S_%f'))).resolve()
    output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(root/'gui'))
    from contour_fields import FieldResults,COMPONENTS
    pack=root/'examples/standard';spec=json.loads((pack/'cases.json').read_text(encoding='utf-8'))
    report={'date_local':now.isoformat(),'timezone':'Asia/Shanghai','python':sys.version,'platform':platform.platform(),'project':'.','build':str(build.relative_to(root)) if build.is_relative_to(root) else str(build),'cases':[]}
    suffix='.exe' if sys.platform=='win32' else ''
    solver=build/('NonlinearFEM'+suffix)
    report['solver_sha256']=hashlib.sha256(solver.read_bytes()).hexdigest() if solver.is_file() else None
    tracked=list((root/'src').glob('*.cpp'))+list((root/'include/nonlinearfem').glob('*.h'))+list((root/'gui').glob('*.py'))+list((root/'tests').glob('*.cpp'))+list((root/'tests').glob('*.py'))+list(pack.glob('*'))+[root/'CMakeLists.txt',Path(__file__).resolve()]
    report['source_sha256']={str(p.relative_to(root)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(tracked) if p.is_file()}
    def execute(command,cwd,log):
        result=subprocess.run(list(map(str,command)),cwd=cwd,capture_output=True,text=True,errors='replace',timeout=180)
        log.write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
        if result.returncode:raise RuntimeError(f'{command[0]} returned {result.returncode}; see {log.name}')
    def read(path):
        with path.open(encoding='utf-8',newline='') as stream:return list(csv.DictReader(stream))
    for case in spec['analytic_cases']:
        entry={'id':case['id'],'status':'FAIL'};folder=output/case['id'];folder.mkdir()
        try:
            text=(pack/case['input']).read_text(encoding='utf-8');(folder/'input.nfem').write_text(text,encoding='utf-8')
            command=[solver]+(['--h8'] if case['dimension']==3 else [])+['input.nfem','solution']
            execute(command,folder,folder/'solver.log');solution=folder/'solution';dim=case['dimension']
            rows=read(solution/'nodes.csv');nodes=[tuple(float(r[k]) for k in ('x','y','z')[:dim]) for r in rows]
            if dim==2:
                actual=[tuple(float(r[k]) for k in ('ux','uy')) for r in rows];reactions=rows
                reaction_keys=('reaction_x','reaction_y')
            else:
                actual=[tuple(float(r[k]) for k in ('ux','uy','uz')) for r in read(solution/'frames.csv')]
                reactions=read(solution/'reactions.csv');reaction_keys=('rx','ry','rz')
            displacement_error=max(abs(a-b) for u,v in zip(actual,case['expected_displacement']) for a,b in zip(u,v))
            if len(actual)!=len(case['expected_displacement']) or not math.isfinite(displacement_error):raise ValueError('Invalid result size or displacement')
            fields=FieldResults.from_folder(solution,nodes,__import__('contour_fields').model_info(solution)[1],[actual])
            indices=(0,1,3) if dim==2 else range(6);field_error=0
            scales={'expected_strain':max(1,max(abs(v) for v in case['expected_strain'])),
                    'expected_PK2':max(case['parameters']['young'],max(abs(v) for v in case['expected_PK2']))}
            for row in fields.recover(0):
                for i in indices:
                    for offset,key in ((0,'expected_strain'),(6,'expected_PK2')):
                        expected=case[key][i];field_error=max(field_error,abs(row[offset+i]-expected)/scales[key])
            force=case['parameters']['force'] if case['kind']=='tension' else 0
            balance=max(abs(sum(float(r[k]) for r in reactions)+(force if d==0 else 0)) for d,k in enumerate(reaction_keys))
            if not all(math.isfinite(float(r[k])) for r in reactions for k in reaction_keys):raise ValueError('Non-finite reactions')
            peak_reaction=max(abs(float(r[k])) for r in reactions for k in reaction_keys)
            entry.update(max_displacement_abs_error=displacement_error,max_field_scaled_error=field_error,reaction_balance_abs_error=balance,peak_reaction=peak_reaction)
            fields.export(folder/'element_fields.csv',0)
            if displacement_error>case['displacement_abs_tol'] or field_error>case['field_scaled_tol'] or balance>max(1,abs(force))*1e-7:raise ValueError('Reference tolerance exceeded')
            if case['kind']=='rotation' and peak_reaction>case['reaction_abs_tol']:raise ValueError('Rigid motion produced reaction')
            entry['status']='PASS'
        except Exception as error:entry['error']=str(error)
        report['cases'].append(entry);print(entry['id'],entry['status'],flush=True)
    for case in spec['existing_case_groups']:
        entry={'id':case['id'],'status':'FAIL'};folder=output/case['id'];folder.mkdir()
        try:
            entry['executable_sha256']={name:hashlib.sha256((build/(name+suffix)).read_bytes()).hexdigest() for name in case['executables']}
            for name in case['executables']:execute([build/(name+suffix)],folder,folder/(name+'.log'))
            entry['status']='PASS';entry['executables']=case['executables']
        except Exception as error:entry['error']=str(error)
        report['cases'].append(entry);print(entry['id'],entry['status'],flush=True)
    report['passed']=all(c['status']=='PASS' for c in report['cases'])
    (output/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# 标准算例实测报告','',f"本地时间：{report['date_local']}",'',f"求解程序 SHA256：`{report['solver_sha256']}`",'', '| 算例 | 状态 | 位移绝对误差 | 场量归一化误差 | 反力平衡误差 |','| --- | --- | --- | --- | --- |']
    for c in report['cases']:lines.append('| '+c['id']+' | '+c['status']+' | '+' | '.join(f"{c[k]:.6g}" if k in c else '见日志' for k in ('max_displacement_abs_error','max_field_scaled_error','reaction_balance_abs_error'))+' |')
    lines+=['','误差定义与通过条件见算例手册；测试程序日志保存在各算例目录。源文件指纹详见 report.json。本报告只对应记录的源码及可执行文件，不自动代表今后的版本。']
    (output/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('Report:',output/'REPORT.md');return 0 if report['passed'] else 1

if __name__=='__main__':sys.exit(main())
