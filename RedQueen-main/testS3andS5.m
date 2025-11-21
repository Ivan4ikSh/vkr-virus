% Testing S26 in S3 Text in Red Queen
global yesgenealogy

yesgenealogy=0;
L=50;
tf=200;
N=10000;
Ub=0.003;
sigmas=[0.05 0.1 0.2]; %average s
% distr='const';
%distr='halfgaussian';
%distr='exponential';
distr='uniform';
beta=10; % sharpness parameter in Good et al for g(s)=exp[(-s/s0)^beta]/s0/gamma(1+1/beta)
% where s0=2*sigma*Gamma(1+1/beta)/Gamma(1+2/beta) to have mean(s)=sigma
runs=3;
maxw2=0; maxadapt=0; maxsast=0;

for ksigma=1:length(sigmas)
    sigma=sigmas(ksigma);
    V=zeros(1,runs); w2=V; ada=V; s_tilde=V;
    for run=1:runs
        [adapt, V_av, Vark_av]=recomb_train2(distr,0,0,sigma,L,N,tf,0,Ub,run);
        V(run)=V_av;            % substitution rate
        ada(run)=adapt;
        w2(run)=Vark_av;        % variance between genomes in allele number
        s_ast(run)=adapt/V_av; % effective selection coefficient/sigma
    end
    
    switch distr
        case 'halfgaussian'
    % Iterating kappa=x0-1
            mm=0;
            kappa=1.5; kappaold=77; stil=1.5; % initial values
            while abs(kappaold-kappa) > 1e-4 && mm < 100
                mm=mm+1;
                kappaold=kappa;
                A=stil*sigma/Ub*sqrt(kappa/(kappa+1));
                stil=sqrt(pi*kappa/(2*kappa+1)*log(A));
                a=log(2*N*Ub*stil*(1+kappa)^1.5/pi/kappa)/log(A);
                kappa=a+sqrt(a^2+a); % corrected 22.03.25
            end
            figure(4)
            subplot(2,2,1)
            plot(pi*kappa/2,mean(ada)/sigma^2,'o')
            text(pi*kappa/2,mean(ada)/sigma^2,sprintf('\sigma=%g',sigma))
            hold on
    
            subplot(2,2,2)
 % corrected 22.03.25
            plot(stil,mean(s_ast)/sigma,'o')
            %text(stil,mean(s_ast)/sigma,sprintf('%g %g',sigma,Ub))
            hold on
        case 'uniform'
           % from Good et al 2012, eq. 21
            coeff=2*gamma(1+1/beta)/gamma(1+2/beta); % s0/sigma to have mean(s)=sigma
            s0=coeff*sigma;
            X=sqrt(beta)*s0/Ub;
            vtosigma2=2*coeff^2*log(N*s0*sqrt(log(N*s0)))/log(X*log(X)^0.5)^2/(1+log(X)/4/log(N*s0))^2;
            figure(4)
            subplot(2,2,1)
            plot(vtosigma2,mean(ada)/sigma^2,'o')
            text(vtosigma2,mean(ada)/sigma^2,sprintf('%g',sigma))
            hold on
            
            subplot(2,2,2)
            % from RR2018 eq. 14 obtained from Good et al eq. 20

            RHS=coeff*(sqrt(2)/beta*log(s0/Ub))^(1/(beta-1)); 
            plot(RHS,mean(s_ast)/sigma,'o')
            text(RHS,mean(s_ast)/sigma,sprintf('%g',sigma))
            hold on
    end
            
    subplot(2,2,3)
    plot(mean(V)/mean(s_ast),mean(w2),'o')
    hold on
    % diagonals
    maxadapt=max(maxadapt,mean(ada)/sigma^2);
    maxsast=max(maxsast,mean(s_ast)/sigma);
    maxw2=max(maxw2,mean(w2));
end % loop in sigma

figure(4)
subplot(2,2,1)
plot([0 1.5*maxadapt],[0 1.5*maxadapt])
% ylabel('V/\sigma')
ylabel('Adaption rate/\sigma ^2')
xlabel('Theory Good et al')
title(sprintf('%s \beta =%g N=%g \n L=%g U_b=%g t_f=%g runs=%g',distr,beta,N,L,Ub,tf,runs))
box off
% hold off
subplot(2,2,2)
plot([0 1.5*maxsast],[0 1.5*maxsast])
ylabel('s*/\sigma')
xlabel('Theory Good et al 2012')
box off
% hold off
subplot(2,2,3)
plot([0 1.5*maxw2],[0 1.5*maxw2])
ylabel('w^2')
xlabel('V/s*')
box off
% hold off