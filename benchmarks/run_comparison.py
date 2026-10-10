"""I compare fixed-lambda B1 workloads, not notebook input/plotting overhead."""
import argparse,csv,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python'))
import numpy as np
import mixed_volume as core

def run(args):
    out=ROOT/'results';out.mkdir(exist_ok=True)
    compiler=shutil.which(args.compiler)
    if compiler is None:raise RuntimeError('C++ compiler not found: '+args.compiler)
    binary=out/'mixed_volume_benchmark'
    command=[compiler,'-std=c++17','-O2',str(ROOT/'cpp/mixed_volume_experiment.cpp'),'-o',str(binary)]
    subprocess.run(command,check=True)
    info=dict(measured_at=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),platform=platform.platform(),
              machine=platform.machine(),python=sys.version,numpy=np.__version__,
              compiler=subprocess.check_output([compiler,'--version'],text=True),compile_command=command,
              settings=vars(args),purpose='Matched B1 workloads; experimental pilot, no error or mixing certificate')
    (out/'benchmark_environment.json').write_text(json.dumps(info,indent=2))
    configs=sorted(set([(args.dimension,l) for l in args.lambdas]+[(n,args.lambda_fixed) for n in args.dimensions]))
    path=out/'comparison.csv'
    fields=['implementation','n','lambda','seed','samples','burn_in','spacing','tolerance','phases','volume_estimate','exact_volume','relative_error','seconds']
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for n,lam in configs:
            for repeat in range(args.repeats):
                seed=args.seed+repeat;exact=core.exact_cube_reference(n,lam)['volume']
                # Alternate ordering to avoid always favoring the second implementation.
                order=['Python','C++'] if repeat%2==0 else ['C++','Python']
                for language in order:
                    print(f'{language}: n={n}, lambda={lam}, repeat={repeat+1}/{args.repeats}',flush=True)
                    if language=='Python':
                        core.TOLERANCE=args.tolerance
                        started=time.perf_counter()
                        result=core.Oracle_B1(n,lam,args.samples,seed=seed,burn_in=args.burn_in,spacing=args.spacing,verbose=False)
                        seconds=time.perf_counter()-started;estimate=result['volume'];phases=len(result['alphas'])
                    else:
                        raw=subprocess.check_output([str(binary),'--benchmark',str(n),str(lam),str(args.samples),str(args.burn_in),str(args.spacing),str(seed),str(args.tolerance)],text=True)
                        estimate,seconds,phases=map(float,raw.strip().split(','))
                    writer.writerow(dict(implementation=language,n=n,**{'lambda':lam},seed=seed,samples=args.samples,burn_in=args.burn_in,spacing=args.spacing,tolerance=args.tolerance,phases=int(phases),volume_estimate=estimate,exact_volume=exact,relative_error=abs(estimate/exact-1),seconds=seconds));f.flush()
    print('Saved',path)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compiler',default='clang++')
    p.add_argument('--dimension',type=int,default=9)
    p.add_argument('--lambdas',type=float,nargs='+',default=[0.5,1.5,3.0])
    p.add_argument('--dimensions',type=int,nargs='+',default=[7,9,11])
    p.add_argument('--lambda-fixed',type=float,default=1.5)
    p.add_argument('--samples',type=int,default=2000)
    p.add_argument('--burn-in',type=int,default=1000)
    p.add_argument('--spacing',type=int,default=3)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--seed',type=int,default=2026)
    p.add_argument('--tolerance',type=float,default=1e-9)
    a=p.parse_args()
    if min(a.dimension,*a.dimensions,a.samples,a.spacing,a.repeats)<1 or a.burn_in<0 or a.seed<0 or not 0<a.tolerance<1 or min(*a.lambdas,a.lambda_fixed)<0:p.error('Invalid benchmark settings')
    run(a)
