"""Experimental mixed-volume oracles for [-1,1]^n + lambda B_2^n.
Sourav Roy's cube-and-ball research experiment.
Implementation assembled in October 2026; this is a numerical experiment related to my earlier research.

Requirements: NumPy. Matplotlib is optional for charts.
Terminal / VS Code: python mixed_volume_experiment.py
Jupyter: %run /full/path/to/mixed_volume_experiment.py

This is an experimental prototype, not a certified ideal-oracle implementation.
Finite hit-and-run burn-in does not certify uniformity. Recorded samples are
correlated. Floating-point chord and label decisions are approximate.
No guaranteed accuracy, mixing time, or polynomial-time FPRAS is claimed.
The code does not assess the thesis proofs.
"""
import math
import time
import numpy as np

# Internal function defaults; the interactive run passes your chosen settings.
M = 1000
SAMPLES = 2000
BURN_IN = 1000
SPACING = 3
THRESHOLD = 0.05
SEED = 2026
TOLERANCE = 1e-9

def validate(n, lam, samples=1, burn_in=0, spacing=1):
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("n must be a positive integer")
    if not np.isfinite(lam) or lam < 0:
        raise ValueError("lambda must be finite and nonnegative")
    for name, value, minimum in [("samples", samples, 1), ("burn_in", burn_in, 0), ("spacing", spacing, 1)]:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")

def projection_cube(x):
    """Euclidean projection onto [-1,1]^n; O(n) arithmetic."""
    return np.clip(np.asarray(x, dtype=float), -1.0, 1.0)

def inside_rounded_cube(x, lam):
    """Numerical membership, using the explicit cube projection."""
    return np.linalg.norm(x - projection_cube(x)) <= lam

def component_label(x):
    """Algorithm 1 classification for the cube; boundaries use float comparisons."""
    return int(np.count_nonzero(np.abs(x) > 1.0))

def hit_and_run_step(x, lam, box_radius, rng, tolerance=None):
    if tolerance is None:
        tolerance = TOLERANCE
    """One approximate chord-uniform transition; not an independent uniform sample."""
    d = rng.normal(size=x.size)
    norm = np.linalg.norm(d)
    while norm == 0:
        d = rng.normal(size=x.size)
        norm = np.linalg.norm(d)
    d /= norm
    mask = d != 0
    t1 = (-box_radius - x[mask]) / d[mask]
    t2 = (box_radius - x[mask]) / d[mask]
    lo = float(np.max(np.minimum(t1, t2)))
    hi = float(np.min(np.maximum(t1, t2)))
    def endpoint(t):
        if inside_rounded_cube(x + t*d, lam):
            return t
        good, bad = 0.0, t
        for _ in range(100):
            if abs(bad-good) <= tolerance:
                break
            mid = (good+bad)/2
            if inside_rounded_cube(x + mid*d, lam):
                good = mid
            else:
                bad = mid
        return good
    lo, hi = endpoint(lo), endpoint(hi)
    return x + rng.uniform(lo, hi)*d

def sample_chain(n, lam, box_radius, samples, rng, burn_in=BURN_IN, spacing=SPACING):
    """Start a new chain at zero and record samples after experimental burn-in."""
    validate(n, lam, samples, burn_in, spacing)
    x = np.zeros(n)
    for _ in range(burn_in):
        x = hit_and_run_step(x, lam, box_radius, rng)
    points = np.empty((samples,n))
    for j in range(samples):
        for _ in range(spacing):
            x = hit_and_run_step(x, lam, box_radius, rng)
        points[j] = x
    return points

class CubeB2:
    """Stateful experimental implementation of Algorithm 1 / B2."""
    def __init__(self, n, lam, seed=SEED, burn_in=BURN_IN, spacing=SPACING):
        validate(n, lam, 1, burn_in, spacing)
        self.n, self.lam, self.spacing = n, lam, spacing
        self.rng = np.random.default_rng(seed)
        self.x = np.zeros(n)
        for _ in range(burn_in):
            self.x = hit_and_run_step(self.x, lam, 1+lam, self.rng)
    def __call__(self):
        for _ in range(self.spacing):
            self.x = hit_and_run_step(self.x, self.lam, 1+self.lam, self.rng)
        return component_label(self.x)

def Oracle_B2(oracle):
    """Return just I in {0,...,n}; the supplied object retains its walk state."""
    return oracle()

def Oracle_B3(n, lam, m=M, seed=SEED, burn_in=BURN_IN, spacing=SPACING):
    """Experimental empirical distribution, using a continuing B2 chain."""
    validate(n, lam, m, burn_in, spacing)
    oracle = CubeB2(n, lam, seed, burn_in, spacing)
    indices = np.empty(m, dtype=int)
    points = np.empty((m,n))
    for j in range(m):
        indices[j] = Oracle_B2(oracle)
        points[j] = oracle.x
    counts = np.bincount(indices, minlength=n+1)
    return dict(probabilities=counts/m, counts=counts, indices=indices, points=points)

def Oracle_B1(n, lam, samples=SAMPLES, phases=None, seed=SEED+10000,
              burn_in=BURN_IN, spacing=SPACING, verbose=True):
    validate(n, lam, samples, burn_in, spacing)
    s0, sf = max(1.0,(1+lam)/math.sqrt(n)), 1+lam
    if s0 == sf:
        return dict(volume=(2*sf)**n, alphas=np.array([]), radii=np.array([sf]))
    if phases is None:
        phases = max(1, math.ceil(n*math.log(sf/s0)/math.log(2)))
    if not isinstance(phases, (int,np.integer)) or phases < 1:
        raise ValueError("phases must be a positive integer")
    radii = np.exp(np.linspace(math.log(s0),math.log(sf),phases+1))
    rng = np.random.default_rng(seed)
    alphas = []
    for j in range(1,phases+1):
        if verbose: print(f"B1 phase {j}/{phases}: sampling...", flush=True)
        points = sample_chain(n,lam,radii[j],samples,rng,burn_in,spacing)
        hits = np.count_nonzero(np.max(np.abs(points),axis=1) <= radii[j-1])
        alpha = hits/samples
        if verbose: print(f"  alpha_{j}={alpha:.6f} ({hits}/{samples})", flush=True)
        if alpha == 0: raise RuntimeError("Zero phase hits; increase samples or phases.")
        alphas.append(alpha)
    alphas = np.asarray(alphas)
    volume = math.exp(n*math.log(2*s0)-np.log(alphas).sum())
    return dict(volume=volume, alphas=alphas, radii=radii)

def exact_mixed_volume(n,k):
    return math.exp((n-k)*math.log(2)+(k/2)*math.log(math.pi)-math.lgamma(k/2+1))

def exact_cube_reference(n,lam):
    validate(n,lam)
    terms = np.array([math.comb(n,i)*exact_mixed_volume(n,i)*lam**i for i in range(n+1)])
    volume = terms.sum()
    return dict(volume=volume, probabilities=terms/volume,
                mixed_volumes=np.array([exact_mixed_volume(n,i) for i in range(n+1)]))

def Algorithm_2(n,k,threshold=THRESHOLD,samples=SAMPLES,seed=SEED,
                max_steps=25,max_refinements=3,burn_in=BURN_IN,spacing=SPACING):
    validate(n,1.0,samples,burn_in,spacing)
    if not isinstance(k,(int,np.integer)) or not 0 <= k <= n: raise ValueError("Invalid k")
    if not 0 < threshold < 1: raise ValueError("threshold must be in (0,1)")
    left, right = math.log(.01/n), math.log(100*n*math.sqrt(n))
    margin, history = threshold/4, []
    for step in range(max_steps):
        mid = (left+right)/2
        lam, decision = math.exp(mid), None
        for refinement in range(max_refinements):
            budget = samples*2**refinement
            p = Oracle_B3(n,lam,budget,seed+100*step+refinement,
                          burn_in*2**refinement,spacing)['probabilities']
            history.append(dict(step=step+1,lam=lam,samples=budget,p=p.copy()))
            print(f"Search {step+1}: lambda={lam:.6g}, samples={budget}, p_{k}={p[k]:.6f}",flush=True)
            if p[k] >= threshold+margin:
                return dict(status='empirical_acceptance',lam=lam,probabilities=p,history=history)
            significant = np.flatnonzero(p >= threshold+margin)
            if significant.size:
                if k < significant[0] and np.all(p[:k+1] < threshold-margin):
                    decision = 'decrease'; break
                if k > significant[-1] and np.all(p[k:] < threshold-margin):
                    decision = 'increase'; break
        if decision == 'decrease': right = mid
        elif decision == 'increase': left = mid
        else: return dict(status='inconclusive_sampling',lam=lam,probabilities=p,history=history)
    return dict(status='search_limit_reached',lam=lam,probabilities=p,history=history)

def Algorithm_3(n,k,samples=SAMPLES,threshold=THRESHOLD,seed=SEED,
                burn_in=BURN_IN,spacing=SPACING):
    validate(n,1.0,samples,burn_in,spacing)
    if not isinstance(k,(int,np.integer)) or not 0 <= k <= n: raise ValueError("Invalid k")
    if k in (0,n):
        return dict(status='exact_cube_endpoint',mixed_volume=exact_mixed_volume(n,k),lam=None)
    search = Algorithm_2(n,k,threshold,samples,seed,burn_in=burn_in,spacing=spacing)
    if search['status'] != 'empirical_acceptance':
        return dict(status=search['status'],mixed_volume=None,search=search)
    lam = search['lam']
    print("Estimating volume...",flush=True)
    b1 = Oracle_B1(n,lam,samples,seed=seed+10000,burn_in=burn_in,spacing=spacing)
    print("Fresh component sampling...",flush=True)
    b3 = Oracle_B3(n,lam,samples,seed+20000,burn_in,spacing)
    pk = b3['probabilities'][k]
    if pk == 0: return dict(status='zero_target_hits',mixed_volume=None,search=search)
    estimate = math.exp(math.log(pk)+math.log(b1['volume'])
                        -math.log(math.comb(n,k))-k*math.log(lam))
    return dict(status='experimental_estimate',mixed_volume=estimate,lam=lam,
                volume=b1['volume'],alphas=b1['alphas'],b3=b3,search=search)

