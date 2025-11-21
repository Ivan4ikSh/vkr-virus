function [adapt V_av Vark_av]=recomb_train2(distribution_s,r,M,s0,L,N,tf,f0,muL,run)
% one can add also TMRCA 
global yesgenealogy ts f1site

%% Arguments:
% distribution_s is one of the four forms below
% r recombination rate per genome
% M crossover number
% s0 average selection coefficient
% L number of loci
% N population number
% tf full time in generations
% f0 initial beneficial allele frequency
% muL mutation rate per genome
% run the number of the run

%% random seed, so can average over runs.
rng(run)

%% Mutation rate per locus
mu=muL/L;                      % time interval for a kind of initial conditions
T= 0:tf;                       % times

%% Distributions of s and parameter titles of fig 1 

switch distribution_s
    case 'exponential'
        %  g(s) = (1/s0)*exp(-s/s0)
        s = -s0*log(rand(1,L));
        ts = sprintf('exp(-s/s_0)\n N=%g,r=%g,L=%g,s_0=%g\n f0=%g, muL=%g, tf=%g, M=%g, run=%g',...
            N,r,L,s0,f0, muL,tf, M, run);
    case 'const'
        s = s0*ones(1,L);
         ts = sprintf('same s=s_0\n N=%g,r=%g,L=%g,s_0=%g\n f0=%g, muL=%g, tf=%g, M=%g, run=%g',...
            N,r,L,s0,f0, muL,tf, M, run);
    case 'halfgaussian'
        % g(s) = (2/pi/s0)*exp(-s^2/pi/s0^2), s > 0, otherwise 
        % mean(s)= s0. If any sign s is allowed, it has std(s) = s0*sqrt(pi/2) 
        s = s0*sqrt(pi/2)*abs(randn(1,L));  % mean(s)=s0 checked directly
        ts = sprintf('halfgaussian\n N=%g,r=%g,L=%g,s_0=%g\n f0=%g, muL=%g, tf=%g, M=%g, run=%g',...
            N,r,L,s0,f0, muL,tf, M, run);
    case 'uniform'
        s=2*s0*rand(1,L);
        ts = sprintf('uniform [0 2 s_0]\n N=%g,r=%g,L=%g,s_0=%g\n f0=%g, muL=%g, tf=%g, M=%g, run=%g',...
            N,r,L,s0,f0, muL,tf, M, run); 
end


%% Initial settings

A=(1:N)'*ones(1,L); % ancestor matrix
Knew=zeros(N,L);  % binary DNA sequences
Anew=zeros(N,L);  % ancestor label sequences
Pnew=zeros(N,L);  % parent labels
%P1=zeros(N,length(T)); Pm=P1; PL=P1; % Initial parent labels for 3 sites in time, for phylogeny 

tint=round(tf/7); % plot time interval 
col='rgbmkrgbmkrgbmkrgbmk';
fsample=0.1; % percent sample for pairs

%% Initial population

% Case 1: randomly distributed good alleles with fixed frequency f0
if f0~=0
    K=(rand(N,L) < f0); 
    % Matrix K: Each row is a genome, sequence of 0 (deleterious) and 1
    % (beneificial)
else
    K=zeros(N,L);
end

P1=zeros(N,tf+1); PL=P1; % Initial parent labels for 3 sites in time, for phylogeny
W=zeros(N,tf+1); % Initial for fitness\

fsite=zeros(tf+1,L); kav=zeros(1,tf+1); Vark=kav; fnot0=kav;
C=kav; Call=kav; meanW=kav;

%% Evolution starts...
for t=T
% beneficial mutation 
    if mu > 0
        K = xor(K, rand(N,L) < mu);
    end
   
    % Parent label initial
    P=(1:N)'*ones(1,L); % parent matrix
%% Recombination of randomly chosen pairs with one-parent replacement
    npairs=round(r*N/2);
    ii=ceil(rand(npairs,2)*N); i1=ii(:,1); i2=ii(:,2);  % 2 columns of random indices of parents
    for i=1:npairs
        % generating random 1 or 0 for each site, with probabilities M/L and 1-M/L,
        % respectively, to mark crossovers with 1
        % Even/odd xx shows site segments copied from 1st parent, 2nd, 1st, 2nd etc
        xx=cumsum(rand(1,L) < M/L); 
        %  sites copied from 1st parent are marked as 1
        first=(round(xx/2)==xx/2);     
        % recombinant DNA sequence
        prog1=K(i1(i),:).*first+K(i2(i),:).*(1-first);   
        % recombinant ancestor label
        prog1A = A(i1(i),:).*first + A(i2(i),:).*(1-first);
        % recombinant parent label
        prog1P = P(i1(i),:).*first + P(i2(i),:).*(1-first);      
%% replacing a parent
        if rand > 0.5  
            K(i1(i),:)=prog1; % 1st parent's DNA replaced
            A(i1(i),:)=prog1A; % 1st parent's ancestor label replaced
            P(i1(i),:)=prog1P; % 1st parent's ancestor label replaced
        else
            K(i2(i),:)=prog1; % 2nd parent's DNA replaced
            A(i2(i),:)=prog1A; % 2nd parent's ancestor label replaced
            P(i2(i),:)=prog1P; % 2nd parent's ancestor label replaced
        end
    end % of mating  
    
%% Random sampling of progeny and natural selection with a broken stick
% column of N log-fitnesses ; state with 0s only has w=0 (fitness 1) by definition
    w=K*(s'); 
% add epistasis here! Input parameters Tij, Eij are set separately in the beginning.
    nprogav=exp(w)/mean(exp(w));    % average progeny number
    b2=cumsum(nprogav);
    b1=[0;b2(1:end-1)];             % broken stick
    X=rand(N,1)*N;
    nprog=zeros(N,1);
    for i=1:N
        nprog(i)=sum( X > b1(i) & X < b2(i)); % actual progeny number
    end
    %disp(sprintf('size(nprog)=%g ',size(nprog)))
  
%% Updating population
    is=[0;cumsum(nprog(1:(N-1)))];
    for i=1:N
        if nprog(i)
            Knew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*K(i,:); % DNA sequences
            Anew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*A(i,:); % ancestor label
            Pnew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*P(i,:); % parent label
        end
    end
    % randomizing numbers of genomes to break families
    %X=randperm(N);
    K=Knew;%(X,:);
    A=Anew;%(X,:);
    P=Pnew;%(X,:);
     
    % control of population preservation
    sK=size(Knew);sK=sK(1);
    if sK  ~= N, disp('N not conserved'),return;end
   
    % updating sequences done
    
    %% Memorizing basic observables 
   
    fsite(t+1,:)=mean(K);        % 1-site allele frequencies at all sites
    kav(t+1)=L*mean(mean(K));     % average allele number per genome
    %Vark(t+1)=(std(w)/s0)^2;    % only for fixed s   
    Vark(t+1)=(std(sum(K,2)))^2;   % variance of the allele number between genomes
    fnot0(t+1)= mean(~all(K==0)); % fraction of loci not uniformly 0
    % Fraction of loci with a common ancestor
    xx=round(N*fsample);
    i1=ceil(N*rand(1,xx)); i2=ceil(N*rand(1,xx));
    C(t+1)=mean(mean(A(i1,:)==A(i2,:))); % pairs
    Call(t+1)= mean(~std(A));            % all population
    % Fitness of all genomes
    W(:,t+1) = w;
    meanW(t+1)=mean(w);
    
    %% memorize names of parents for 1st and last loci   
    P1(:,t+1)=P(:,1); PL(:,t+1)=P(:,L);    
    
    %% Plotting the wave and ancestral clone  spectrum at some time points
    if tint*round(t/tint)==t 
        c=col(round(t/tint)+1);
        % K(1:6,1:8)
        figure(1)
     subplot(2,2,1)
        [nn,xx]=hist(w);         % histogram of fitness among genomes
        plot(xx/s0,nn,c)
        hold on
    end % if sampling time 
end % loop in time

% Evolution has ended

%% Final plots

figure(1) % Traveling wave
subplot(2,2,1)
hold off
xlabel('Allele number,  k  ');
ylabel('Density distribution')
title(ts)   % title  with parameter values
 
%% Average observables vs time

kav=reshape(kav, size(T)); % average k
Vark=reshape(Vark, size(T)); % Variance k
dist=mean(L*2*fsite.*(1-fsite),2);
dist=reshape(dist,size(T));
% calculating instantaneous speed averaged over 3 generations
V=(kav(4:end)-kav(1:(end-3)))/3; 
V=[V(1:3),V];
%Fisher = V./Vark/s0; % checking Fisher theorem

% 1-site theory without mutation
if f0 
    t50=(1/s0)*log(f0/(1-f0));
    f1site=1-(1+exp(s0*(T-t50))).^(-1); 
else
    t50=(1/s0)*log(s0/mu);
    f1site=1-(mu/s0+(1+exp(s0*(T-t50))).^(-1)); 
end


%%
subplot(2,2,2)
plot(T,kav/L,'r',T,sqrt(Vark)/L,'b',T,dist/L,'g',T,C,'m',...
    T,fnot0,'k',T,Call,'m--')%,T,Fisher,'oc');
title(sprintf('k_{av}/L r SD_k/L b, dist/L g \n C m, f_{not0} k, Call m- -'))% Fisher c'))
xlabel('Time, t');
axis([0 T(end) 0 1.2])
Cinf=fnot0(end);

subplot(2,2,3)
avfend=mean(fsite(end,:));
% plot(T',fsite(:,[1,10:10:L]),T,f1site,'k--')
plot(T',fsite,T,f1site,'k--')
ylabel('Allele frequency at all sites')
xlabel ('Time, t')
title(sprintf('avfend = %g',avfend))

subplot(2,2,4)
for t=1:tint:tf
    [nn,xx]=hist(fsite(t,:));
    plot(xx,nn)
    hold on
end
ylabel('Histogram at a few fixed times')
xlabel ('Allele frequency at a site, t')
hold off


%% Adaptation rate averaged over the last 3/4 of time interval
t1=round(tf/4+1); 
adapt=(meanW(tf+1) - meanW(t1))/(tf-t1);

%% Substitution rate averaged over the last 3/4 of time interval
V_av=(kav(tf+1) - kav(t1))/(tf-t1);

%% Variance of allele number averaged over the last 3/4 of time interval
Vark_av=mean(Vark(t1:(tf+1)));

if yesgenealogy
%% Genealogy trajectory
figure(2)
% Sample of genomes at t=tf
ii=[1 N/4:N/4:N];
% Genealogy numbers loci 1 and 2 
G1(tf+1,ii) = ii; 
GL(tf+1,ii) = ii;
for t = tf:-1:1
    G1(t,ii) = P1(G1(t+1,ii),t+1);
    GL(t,ii) = PL(GL(t+1,ii),t+1);
end
subplot(2,2,1)
plot(T,G1(:,ii)); 
title(ts)
ylabel('Parent, locus 1')

%% Calculating TMRCA 
TMRCA= tf - find(~std(G1(:,ii),2),1,'last');
%%

subplot(2,2,3)
plot(T,GL(:,ii)); 
xlabel('Time')
ylabel('Parent, locus L')

%% Fitness trajectory
w1=zeros(tf+1,length(ii)); wL=zeros(tf+1,length(ii));
for t = T
    w1(t+1,ii)=W(G1(t+1,ii),t+1);
    wL(t+1,ii)=W(GL(t+1,ii),t+1);
end
subplot(2,2,2)
plot(T,w1)
xlabel('Time')
ylabel('Fitness')
subplot(2,2,4)
plot(T,wL)
xlabel('Time')
ylabel('Fitness')

%% Phylogenetic tree
figure(3)
m=length(ii);
yesconnect1=ones(size(ii)); % plotting mark
yesconnectL=ones(size(ii)); % plotting mark
color='rbmkrbmkrbmk'; % colors of lineage

% Loop back in time 
for t=tf:-1:0
    % first site tree
     subplot(2,1,1)
    for i=1:m 
        % find all intersections > i
        jj=find(G1(t+1,ii(i))==G1(t+1,ii((i+1):m)));
        if isempty(jj) && yesconnect1(i)
            % if none, plot a horizontal segment
            plot ([t t+1],[i i],color(i)); hold on 
        elseif yesconnect1(i)
            % plot a segment angled up
            plot ([t t+1],[max(jj)+i, i],color(i)); hold on 
            % and do not plot it anymore
            yesconnect1(i)=0;  
        end 
    end % in lineages
    
    % last site tree
    subplot(2,1,2)
    for i=1:m 
        % find all intersections > i
        jj=find(GL(t+1,ii(i))==GL(t+1,ii((i+1):m)));
        if isempty(jj) && yesconnectL(i)
            % if none, plot a horizontal segment
            plot ([t t+1],[i i],color(i)); hold on 
        elseif yesconnectL(i)
            % plot a segment angled up
            plot ([t t+1],[max(jj)+i, i],color(i)); hold on 
            % and do not plot it anymore
            yesconnectL(i)=0;  
        end 
    end % in lineages
end  % in time

subplot(2,1,1)
hold off
ylabel('Tree, locus 1')
xlabel('Time')
title(ts)
axis([0 tf+1 0.5 m+0.5])
subplot(2,1,2)
hold off
ylabel('Tree, locus L')
xlabel('Time')
title(ts)
axis([0 tf+1 0.5 m+0.5])
end % yesgenealogy?
% end program
