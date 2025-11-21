%Made with MATLAB2022a
%% parameters
% generate (inside the s_measure program - gen) - if you want to test the programm set it to true
% parameters required to work with Monte Carlo simulation:
% run - number of Monte Carlo evolution runs, also used as the seed for the
% first run
% r - recombination probability per genome
% M - crossover number
% L - number of loci
% N - population size
% run2 - the seed for generating the random distribution of s
% NUbs2 - N*Ub*s^2 (Ub - probability of the benefitial mutation)
% tf - end time of evolution
% s0 - the width of the uniform distribution of selection coefficient
% f0 - initial value of f
% ac*s0, bc*s0 - borders of uniform s distribution
% mu - mutation probability per site
% tsec1, tsec2, tsec3 - time of the first, second and third sequence sample
% parameters required to work both with real and generated data
% C - initial value of C in equation 1
% appr - method of curve approximation:
% 'poly' - by basic polynomials
% 'spline' - by cubic splines
% 'pchip' - by Piecewise Cubi Hermite Polynomial (PCHIP)
% 'test' - without any approximation and C finding (only for curves printing)
% apprR - additional parameter for 'poly' approximation - rate of
% polynomial

global tsec1 tsec2 tsec3 r s0 M L N tf f0 run run2 mu ac bc NUbs2


generate=false;

run=100;
r=0;
M=2;
L=100;
N=1000;
run2=1;
NUbs2=0.25;
tf=60;
    s0=0.05;
    f0=0.1;
    ac=-1;
    bc=1;
    mu=NUbs2/N/s0^2/L;
    tsec1=20;
    tsec2=40;
    tsec3=60;
    C=0;
    appr='test'; 
    apprR=4;

if generate
    
for i =1:run
[dat{i},~,~,~,~,~,~]=recomb_2022(r,s0,ac,bc,M,L,N,tf,f0,i,run2);
gen1{i}=dat{i}{1};
gen2{i}=dat{i}{2};
gen3{i}=dat{i}{3};
end
    
order=1:L;
sdis = s_measure(gen1.',gen2.',gen3.',order,C,appr,apprR,generate);
else
    load('BinData.mat','data');
genome1 = data{1,1};
genome2 = data{1,2};
genome3 = data{1,3};
order = data{1,4};
   sdis = s_measure(genome1,genome2,genome3,order,C,appr,apprR,generate);
end


