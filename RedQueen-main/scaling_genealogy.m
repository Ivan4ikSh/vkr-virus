% scaling of TMRCA
% TMRCA=z*sqrt(2*log(N*s)/v);
global yesgenealogy

yesgenealogy=1;
L=100;
tf=400;
N=10000;
Ub=0.001;
% distr='const';
% distr='halfgaussian';
distr='exponential';
runs=20;
TT=zeros(1,runs); ada=TT;
for s0=[0.025 0.05 0.1 0.2]
    for run=1:runs
        [TMRCA, adapt]=recomb_train(distr,0,0,s0,L,N,tf,0,Ub,run)
        if ~isempty(TMRCA)
            TT(run)=TMRCA; 
        else
            TT(run)=TT(run-1);
        end
        ada(run)=adapt;
    end
    figure(4)
    plot(sqrt(2*log(N*s0)/mean(ada)),mean(TT),'o')
    hold on
end
ylabel('T_{MRCA}')
xlabel('[2log(Ns_0)/v]^{1/2}')
title(sprintf('%s N=%g Ub=%g L=%g tf=%g',distr,N,Ub,L,tf))
