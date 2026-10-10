// Experimental cube + Euclidean ball oracles, not certified ideal oracles.
// I use finite hit-and-run burn-in, correlated samples, and approximate chords.
// I do not claim a mixing-time bound or guaranteed approximation error.
// Compile: clang++ -std=c++17 -O2 mixed_volume_experiment.cpp -o mixed_volume
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
using Vec = std::vector<double>;
struct Settings { int n, k; long long m, samples, burn, spacing; double lam, threshold, tol; unsigned long long seed; };

double ask(const std::string& explanation, const std::string& name, double minimum,
           bool integer=false, bool strict=false, double maximum=1e15) {
    std::cout << "\n" << explanation << "\n";
    for (;;) {
        std::cout << name << ": " << std::flush;
        std::string line; if (!std::getline(std::cin,line)) throw std::runtime_error("Input ended.");
        try {
            size_t used; double v=std::stod(line,&used);
            if (line.find_first_not_of(" \t\r",used)!=std::string::npos || !std::isfinite(v)
                || v<minimum || (strict && v==minimum) || v>maximum
                || (integer && std::floor(v)!=v)) throw std::invalid_argument("value");
            return v;
        } catch (...) { std::cout << "Enter a valid " << (integer?"integer":"number") << " in the stated range.\n"; }
    }
}

// I compute distance to the cube by coordinate clipping.
bool inside(const Vec& x,double lam) {
    double d2=0; for(double v:x) { double d=v-std::clamp(v,-1.0,1.0); d2+=d*d; }
    return std::sqrt(d2)<=lam;
}
int label(const Vec& x) { int i=0; for(double v:x) i+=std::abs(v)>1.0; return i; }

// I draw a sphere direction, bracket the chord in a containing box,
// and retain inside endpoints after numerical bisection.
void step(Vec& x,double lam,double box,std::mt19937_64& rng,double tol) {
    std::normal_distribution<double> normal(0,1); Vec d(x.size()); double norm;
    do { norm=0;for(double& v:d){v=normal(rng);norm+=v*v;}norm=std::sqrt(norm); } while(norm==0);
    for(double& v:d)v/=norm;
    double lo=-std::numeric_limits<double>::infinity(),hi=-lo;
    for(size_t j=0;j<x.size();++j) if(d[j]!=0) {
        double a=(-box-x[j])/d[j],b=(box-x[j])/d[j];
        lo=std::max(lo,std::min(a,b));hi=std::min(hi,std::max(a,b));
    }
    auto member=[&](double t){Vec y(x.size());for(size_t j=0;j<x.size();++j)y[j]=x[j]+t*d[j];return inside(y,lam);};
    auto endpoint=[&](double t){
        if(member(t))return t;
        double good=0,bad=t;
        for(int q=0;q<100 && std::abs(bad-good)>tol;++q){double mid=(good+bad)/2;if(member(mid))good=mid;else bad=mid;}
        return good;
    };
    lo=endpoint(lo);hi=endpoint(hi);
    double t=std::uniform_real_distribution<double>(lo,hi)(rng);
    for(size_t j=0;j<x.size();++j)x[j]+=t*d[j];
}

// Algorithm 1 / B2: I keep the previous point between queries.
class B2 {
    Vec x; double lam,tol;long long spacing; std::mt19937_64 rng;
public:
    B2(const Settings& s,double lambda,unsigned long long seed,long long burn)
      :x(s.n,0),lam(lambda),tol(s.tol),spacing(s.spacing),rng(seed){for(long long j=0;j<burn;++j)step(x,lam,1+lam,rng,tol);}
    int operator()(){for(long long j=0;j<spacing;++j)step(x,lam,1+lam,rng,tol);return label(x);}
};

// B3: I return empirical label probabilities from a continuing B2 chain.
Vec b3(const Settings& s,double lam,long long m,unsigned long long seed,long long burn) {
    B2 oracle(s,lam,seed,burn);Vec p(s.n+1,0);
    for(long long j=0;j<m;++j)++p[oracle()];
    for(double& v:p) { v/=m; }
    return p;
}
double logchoose(int n,int k){return std::lgamma(n+1.0)-std::lgamma(k+1.0)-std::lgamma(n-k+1.0);}
double exactmixed(int n,int k){return std::exp((n-k)*std::log(2.0)+k*.5*std::log(std::acos(-1.0))-std::lgamma(k*.5+1));}
Vec exactterms(int n,double lam){Vec a(n+1);for(int i=0;i<=n;++i)a[i]=(lam==0?(i==0?exactmixed(n,0):0):std::exp(logchoose(n,i)+std::log(exactmixed(n,i))+i*std::log(lam)));return a;}
void printvec(const Vec& v){std::cout<<"[";for(size_t j=0;j<v.size();++j)std::cout<<(j?", ":"")<<v[j];std::cout<<"]\n";}
struct Volume{double value;Vec alphas;};

// B1: I estimate nested-body ratios and multiply their reciprocals.
Volume b1(const Settings& s,double lam,unsigned long long seed) {
    double sf=1+lam,s0=std::max(1.0,sf/std::sqrt(double(s.n)));
    if(s0==sf)return {std::pow(2*sf,s.n),{}};
    int phases=std::max(1,int(std::ceil(s.n*std::log(sf/s0)/std::log(2.0))));
    double logv=s.n*std::log(2*s0);Vec alphas;std::mt19937_64 rng(seed);
    for(int j=1;j<=phases;++j){
        double prev=s0*std::exp(std::log(sf/s0)*(j-1)/phases),box=s0*std::exp(std::log(sf/s0)*j/phases);
        std::cout<<"B1 phase "<<j<<"/"<<phases<<": sampling...\n"<<std::flush;
        Vec x(s.n,0);for(long long t=0;t<s.burn;++t)step(x,lam,box,rng,s.tol);
        long long hits=0;
        for(long long t=0;t<s.samples;++t){
            for(long long q=0;q<s.spacing;++q)step(x,lam,box,rng,s.tol);
            bool hit=true;for(double v:x)if(std::abs(v)>prev){hit=false;break;}hits+=hit;
        }
        double a=double(hits)/s.samples;
        std::cout<<"  alpha_"<<j<<" = "<<a<<" ("<<hits<<"/"<<s.samples<<")\n";
        if(a==0)throw std::runtime_error("Zero phase hits; increase the sampling budget.");
        alphas.push_back(a);logv-=std::log(a);
    }
    return {std::exp(logv),alphas};
}

// Algorithm 2: I use an experimental bracket and threshold margin.
// Inconclusive observations cause a refusal to return a forced estimate.
double algorithm2(const Settings& s){
    double left=std::log(.01/s.n),right=std::log(100*s.n*std::sqrt(double(s.n))),margin=s.threshold/4;
    for(int t=0;t<25;++t){
        double mid=(left+right)/2,lam=std::exp(mid);int decision=0;
        for(int r=0;r<3;++r){
            auto p=b3(s,lam,s.samples*(1LL<<r),s.seed+100*t+r,s.burn*(1LL<<r));
            std::cout<<"Search "<<t+1<<": lambda="<<lam<<", p_"<<s.k<<"="<<p[s.k]<<"\n"<<std::flush;
            if(p[s.k]>=s.threshold+margin)return lam;
            int first=-1,last=-1;for(int i=0;i<=s.n;++i)if(p[i]>=s.threshold+margin){if(first<0)first=i;last=i;}
            if(first>=0 && s.k<first){bool low=true;for(int i=0;i<=s.k;++i)low &= p[i]<s.threshold-margin;if(low){decision=-1;break;}}
            if(last>=0 && s.k>last){bool low=true;for(int i=s.k;i<=s.n;++i)low &= p[i]<s.threshold-margin;if(low){decision=1;break;}}
        }
        if(decision<0)right=mid;else if(decision>0)left=mid;else throw std::runtime_error("Search inconclusive: no mixed-volume estimate returned.");
    }
    throw std::runtime_error("Search limit reached.");
}
// Algorithm 3: I combine B1 with a fresh B3 query at the selected lambda.
double algorithm3(const Settings& s){
    if(s.k==0 || s.k==s.n){std::cout<<"Exact cube endpoint case.\n";return exactmixed(s.n,s.k);}
    double lam=algorithm2(s);auto volume=b1(s,lam,s.seed+10000);
    auto p=b3(s,lam,s.samples,s.seed+20000,s.burn);
    std::cout<<"Selected lambda: "<<lam<<"\nFresh B3 vector: ";printvec(p);
    if(p[s.k]==0)throw std::runtime_error("No target hits in fresh B3 query.");
    return std::exp(std::log(p[s.k])+std::log(volume.value)-logchoose(s.n,s.k)-s.k*std::log(lam));
}
// Machine-readable benchmark: n lambda samples burn spacing seed tolerance.
// The timer includes only B1 and suppresses progress output.
int main(int argc,char** argv){try{
    if(argc==9 && std::string(argv[1])=="--benchmark") {
        Settings s{};
        s.n=std::stoi(argv[2]);s.lam=std::stod(argv[3]);s.samples=std::stoll(argv[4]);
        s.burn=std::stoll(argv[5]);s.spacing=std::stoll(argv[6]);s.seed=std::stoull(argv[7]);s.tol=std::stod(argv[8]);
        if(s.n<1 || s.lam<0 || !std::isfinite(s.lam) || s.samples<1 || s.burn<0 || s.spacing<1 || s.tol<=0 || !std::isfinite(s.tol))
            throw std::runtime_error("Invalid benchmark settings.");
        auto* buffer=std::cout.rdbuf();std::cout.rdbuf(nullptr);
        auto start=std::chrono::steady_clock::now();auto v=b1(s,s.lam,s.seed);
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        std::cout.rdbuf(buffer);
        std::cout<<std::setprecision(17)<<v.value<<","<<seconds<<","<<v.alphas.size()<<"\n";
        return 0;
    }
    if(argc!=1)throw std::runtime_error("Use interactive mode or --benchmark n lambda samples burn spacing seed tolerance.");
    std::cout<<std::setprecision(10)<<"Experimental mixed-volume oracles for a centered cube plus a Euclidean ball.\n";
    Settings s;
    s.n=int(ask("Dimension is the number of coordinates per point.","Choose n (positive integer)",1,true,false,100000));
    s.lam=ask("Lambda is the radius of the ball added to the cube; zero leaves the cube unchanged.","Choose lambda",0);
    s.k=int(ask("The target index labels the mixed-volume coefficient; 0 means cube volume and n means ball volume.","Choose k (0 through n)",0,true,false,s.n));
    s.m=(long long)ask("This is the number of labels recorded for the standalone B3 distribution.","Choose m",1,true);
    s.samples=(long long)ask("This separate budget is used per B1 phase and initial Algorithm 2 query.","Choose samples",1,true);
    s.burn=(long long)ask("Burn-in is the number of initial moves discarded; it is not a certified mixing time.","Choose burn-in",0,true);
    s.spacing=(long long)ask("Spacing is the number of moves before each recorded point; it does not guarantee independence.","Choose spacing",1,true);
    s.threshold=ask("The threshold is a desired component probability for the search, not an error tolerance.","Choose threshold (strictly between 0 and 1)",0,false,true,std::nextafter(1.0,0.0));
    s.seed=(unsigned long long)ask("A seed controls the random sequence; C++ samples need not match Python samples.","Choose nonnegative seed",0,true);
    s.tol=ask("Chord tolerance controls the bisection stopping width; smaller is limited by floating-point precision.","Choose positive tolerance",0,false,true);
    auto start=std::chrono::steady_clock::now();
    B2 oracle(s,s.lam,s.seed,s.burn);std::cout<<"\nB2 label I: "<<oracle()<<"\n";
    auto p=b3(s,s.lam,s.m,s.seed,s.burn);std::cout<<"B3 empirical distribution: ";printvec(p);
    auto terms=exactterms(s.n,s.lam);double exactvolume=0;for(double v:terms)exactvolume+=v;
    for(double& v:terms) { v/=exactvolume; }
    std::cout<<"Exact cube probabilities: ";printvec(terms);
    auto volume=b1(s,s.lam,s.seed+10000);std::cout<<"All alpha values: ";printvec(volume.alphas);
    std::cout<<"Estimated volume: "<<volume.value<<"\nExact cube volume: "<<exactvolume<<"\nObserved relative error: "<<std::abs(volume.value/exactvolume-1)<<"\n";
    try{double estimate=algorithm3(s),exact=exactmixed(s.n,s.k);
        std::cout<<"Estimated M_k: "<<estimate<<"\nExact cube M_k: "<<exact<<"\nObserved relative error: "<<std::abs(estimate/exact-1)<<"\n";
    }catch(const std::exception& e){std::cout<<"Algorithm 3 status: "<<e.what()<<"\n";}
    std::cout<<"Elapsed seconds: "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\nExperimental results; no certified mixing or accuracy guarantee.\n";
    return 0;
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<"\n";return 1;}}
