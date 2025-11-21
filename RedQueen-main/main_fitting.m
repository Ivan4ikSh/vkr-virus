% Fitting b_y and U_b to data, no vaccine

global distribution_s R0 N TMRCA_exp V_exp du py TMRCA z V finf u K r Kder sigma x0 s_ast

%% Basic parameters

% distribution_s = 'constant'; 
% distribution_s ='exponential';
%distribution_s ='halfgaussian';
distribution_s='uniform';

N = 1e8;           % total population size
z=3; 

%% Experimental values
% Influenza A H2N3
R0 = 1.8;          % basic reproduction number without vaccine
TMRCA_exp =  220;   % [cycle] 3.0 year, trec = 5 day per cycle
V_exp = 0.036;      % 1/[cycle] 2.6 aa/year
log10Ub = (-4:0.1:-2);
by = 2:0.05:10;

% % Influenza A H1N1
% R0=1.5; %H1N1
% TMRCA_exp =  336;   % [cycle] 4.6 year
% V_exp = 0.031;      % 1/[cycle] 2.3 aa/year
% log10Ub = (-4:0.1:-2.2);
% by = 4:0.05:12;

py=0.024;       % susceptibility to homologous strain
du=0.2;         % step in u

SD2=zeros(length(log10Ub),length(by));
TMRCAs=SD2; Vs=SD2; finfs=SD2; sigmas=SD2; x0s=SD2; s_asts=SD2;
err=zeros(1,length(log10Ub));
n=zeros(1,length(log10Ub));

% memorizing observables

for i = 1:length(log10Ub)
    for j = 1:length(by)
       SD2(i,j) = call_fitting([log10Ub(i), by(j)]);
       TMRCAs(i,j)=TMRCA;
       Vs(i,j)=V;
       finfs(i,j)=finf;
       sigmas(i,j)=sigma;
       x0s(i,j)=x0;
       s_asts(i,j)=s_ast;
    end
    [err(i), jmin(i)]=min(SD2(i,:));
end

% location of best fit
[errmin, imin] = min(err);
disp([log10Ub(imin), by(jmin(imin)), errmin])

figure(1)
[C,h]=contour(by,log10Ub,real(SD2),[0.005 0.02  0.05 0.1 0.2 0.5 1.0]);
clabel(C,h)
xlabel('b_y')
ylabel('log_{10}U_b')
title(sprintf(['%s R0=%g N=%g py=%g z=%g V_{exp}=%g TMRCA_{exp}=%g \n Best fit: log10Ub=%g  by=%g  SD2=%0.3f'],distribution_s,R0, N,py,z,V_exp,TMRCA_exp,log10Ub(imin), by(jmin(imin)), errmin))
text(by(jmin(imin)),log10Ub(imin),'+')
box off

figure(2)
[C,h]=contour(by,log10Ub,real(TMRCAs/TMRCA_exp),[0.5 0.75 1 1.5 2 3] );
clabel(C,h)
% xlabel('b_y')
% ylabel('log_{10}U_b')
% title('T_{MRCA} relative to experimental value')
% text(by(jmin(imin)),log10Ub(imin),'+')
box off

%figure(3)
hold on
[C,h]=contour(by,log10Ub,real(Vs/V_exp),[0.3 0.5 1 2 3 4 6]);
clabel(C,h)
xlabel('b_y')
ylabel('log_{10}U_b')
title('Speed V and TMRCA relative to observed values')
text(by(jmin(imin)),log10Ub(imin),'+')
box off

figure(4)
[C,h]=contour(by,log10Ub,real((365/5)*finfs), [0.01 0.02 0.05 0.1 0.2 0.4]);
clabel(C,h)
xlabel('b_y')
ylabel('log_{10}U_b')
title('Annual incidence')
text(by(jmin(imin)),log10Ub(imin),'+')
box off

figure(5)
[C,h]=contour(by,log10Ub,sigmas);
clabel(C,h)
xlabel('b_y')
ylabel('log_{10}U_b')
title('Selection coefficient \sigma')
text(by(jmin(imin)),log10Ub(imin),'+')
box off


% figure(6)
% [C,h]=contour(by,log10Ub,real(x0s));
% clabel(C,h)
% xlabel('b_y')
% ylabel('log_{10}U_b')
% title('Lead x_0')
% text(by(jmin(imin)),log10Ub(imin),'+')
% box off


figure(7)
[C,h]=contour(by,log10Ub,real(s_asts));
clabel(C,h)
xlabel('b_y')
ylabel('log_{10}U_b')
title('Effective s')
text(by(jmin(imin)),log10Ub(imin),'+')
box off

figure(6)
[C,h]=contour(by,log10Ub,real(Vs./s_asts),[0.2 0.3 0.5 1 2]);
clabel(C,h)
xlabel('b_y')
ylabel('log_{10}U_b')
title('Peak variance w^2')
text(by(jmin(imin)),log10Ub(imin),'+')
box off

call_fitting([log10Ub(imin), by(jmin(imin))]);
%% Plot K(u)
figure(8)
plot(u,K,u,r,u,Kder)
ylabel('K     Kder     r')
xlabel('Antigenic distance  u')
title(sprintf(['%s R0=%g N=%g py=%g z=%g V_{exp}=%g TMRCA_{exp}=%g \n Best fit: log10Ub=%g  by=%g  SD=%0.3f'],distribution_s,R0, N,py,z,V_exp,TMRCA_exp,log10Ub(imin), by(jmin(imin)), errmin))

figure(1)

%% Automatic search of minimum: very slow
% log10Ub_0 = -4;       % initial approximation for mutation rate per epitope per transmission
% by_0 = 3;        % initial approximation for steepness of K(u) 
% optimset('Display', 'iter','MaxIter',100,'TolX',0.05);
% X = fminsearch(@call_fitting,[log10Ub_0 by_0]);
% 
% Ub=exp(X(1))
% by=X(2)
% disp([TMRCA/TMRCA_exp V/V_exp])


