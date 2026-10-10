"""I plot observed runtime and observed error in two panels for each sweep."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def plot_sweep(rows,kind,fixed,save=True):
    xkey='lambda' if kind=='lambda' else 'n'
    fixedkey='n' if kind=='lambda' else 'lambda'
    selected=[r for r in rows if float(r[fixedkey])==fixed]
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for language,color in [('Python','#2563eb'),('C++','#d97706')]:
        subset=[r for r in selected if r['implementation']==language]
        xs=sorted(set(float(r[xkey]) for r in subset))
        for ax,key,scale in [(axes[0],'seconds',1),(axes[1],'relative_error',100)]:
            vals=[np.array([float(r[key])*scale for r in subset if float(r[xkey])==x]) for x in xs]
            med=np.array([np.median(v) for v in vals]);low=np.array([v.min() for v in vals]);high=np.array([v.max() for v in vals])
            ax.errorbar(xs,med,yerr=[med-low,high-med],marker='o',capsize=4,label=language,color=color)
            ax.set_xlabel('Lambda' if kind=='lambda' else 'Dimension n');ax.grid(alpha=.25);ax.legend()
    axes[0].set_yscale('log');axes[0].set_ylabel('B1 calculation seconds (log scale)');axes[0].set_title('Measured runtime')
    axes[1].set_ylabel('Absolute relative volume error (%)');axes[1].set_title('Observed error against exact cube volume');axes[1].set_ylim(bottom=0)
    fig.suptitle(f'Fixed {fixedkey}={fixed:g}; median and full observed range across seeds')
    fig.tight_layout()
    if save:fig.savefig(ROOT/'results'/f'{kind}_comparison.svg',bbox_inches='tight')
    return fig

def load_rows():
    with (ROOT/'results/comparison.csv').open() as f:return list(csv.DictReader(f))
if __name__=='__main__':
    import json
    settings=json.loads((ROOT/'results/benchmark_environment.json').read_text())['settings']
    rows=load_rows();plot_sweep(rows,'lambda',settings['dimension']);plot_sweep(rows,'dimension',settings['lambda_fixed']);plt.show()
