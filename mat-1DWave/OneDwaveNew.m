% Model of 1D wave with reinfection with ordered variants 
% based on Lin et al JTB 2003, but mutations in discrete space diffusion
% global  R0 a D
% R0 is BRR, D is rescaled mutation rate, a is half-width of crossimmunity 

%% Parameter values without vaccine
 R0=2.5;
 Ub=1e-3; 
 a=7;%14;      % crossimmunity/width for infection 
 N=1e8;    % population size for Poisson tail
 L=200;  xx=(0:L-1)/a;  K=xx./(1+xx);%1-exp(-xx);  % crossimmunity function
 run=1;     % random generator seed 
 rng(run)

%% Vaccine parameters
 %typevac='last';  % the last vaccine only
 typevac='all ';  % all vaccines
 az=4;      % response width for vaccine 
 Dz=2;      % vaccine lag behind virus peak at vaccination
 tv=28;     % vaccination interval
 ez=0.4;    % maximal efficacy
 cz=0;    % vaccine coverage 
%%

% Technical parameters
Tmax=2000;     % Max time. One unit = 1/nu ~ 1 wk
T0=1000;        % Start show, vaccine, turn ordering x < y off
ts=sprintf('K=x/(1+x) %s R0=%g Ub=%g a=%g T0=%g Tmax=%g N=%g \n az=%g Dz=%g tv=%g cz=%g ez=%g run=%g',typevac, R0,Ub,a,T0,Tmax,N,az,Dz,tv,cz,ez,run);
stept=0.5;     % step in time
stepshow=Tmax/5;
M=Tmax/stept;       % number of time points

% Initial conditions
I=zeros(1,L); I(20)=1e-3;             % Infected
R=zeros(1,L); R(1:19)=(1-1e-3)/19;    % Recovered from previous infection, as S
Q=zeros(1,L); P=zeros(1,L);
Inew=zeros(1,L);
figure(1);clf
norm=zeros(1,M);
finf=zeros(1,M);
z=sum((0:L-1).*I)/sum(I);
meanI=zeros(1,M);
Re = R0*ones(1,L);
% Loop in time

for i=1:M
    meanI(i)=sum((0:L-1).*I)/sum(I);
    if stept*i > T0 && round(stept*i/tv)==stept*i/tv
        z = meanI(i) - Dz;
        xx = abs(((0:L-1) - z))/az;
        Kz = 1 - exp(-xx);
        if typevac=='last'
            Re = R0*(1-ez*cz*(1-Kz)); % effective reproduction number w/vaccine
        elseif typevac=='all '
            Re = Re.*(1-ez*cz*(1-Kz));
        end
        figure(1)
        if cz>0
            plot(0:L-1, 0.01*(1-Kz),'g')
            hold on
        end
    end
 
    % Loop in x
    norm(i)=sum(I)+sum(R);  % should be = 1
    finf(i)=sum(I)/(sum(I)+sum(R)); % fraction of infected
    for j=1:L
        Q(j)=sum(R(1:j).*K(j:-1:1));
        P(j)=sum(Re(j:L).*I(j:L).*K(1:(L-j+1)));
        Inew(j)=I(j)*(1+stept*(Re(j)*Q(j)-1)); 
    % genetic drift
        if Inew(j)*N < 10
            xm=round(6*Inew(j)*N);
    % Poisson, broken stick
            Inew(j)=sum(xm*rand(1,xm) < Inew(j)*N)/N;
        end 
        R(j)=R(j)*(1-stept*P(j))+stept*I(j);
    % mutation term
        if j > 1 & j < L
            in=stept*Ub*(I(j+1)+I(j-1));
            out=stept*Ub*2*I(j);
            if in*N < 10
                xm=round(6*in*N);
            % Poisson, broken stick
                in=sum(xm*rand(1,xm) < in*N)/N;
            end 
            Inew(j)=Inew(j)+in-out;
        end
    end
    Inew(1)=I(1)+stept*Ub*(I(2)-I(1));
    Inew(L)=I(L)+stept*Ub*(I(L-1)-I(L));
    I=Inew/norm(i);  
    % plotting
    if round(stept*i/stepshow)==stept*i/stepshow & stept*i > T0
        figure(1)
%         plot(0:L-1,10*I,'r',0:L-1,R,'b')
    semilogy(0:L-1,10*I,'r',0:L-1,R,'b')
    axis([0 L 1e-6 1])
        hold on
    end
    % reference position for speed
    if stept*i==T0
        meanI0=meanI(i);
    end
%     QI(i)=sum(Q.*I);
%     PR(i)=sum(P.*R);
end % loop in time
hold off
xlabel('Mutation number')
ylabel('R (b) 10*I (r) (1-Kz)/100 (m)')
title(ts)

% figure(2)
% % normalization and steady process control
% plot(stept*(1:M),norm,'r',stept*(1:M),meanI,'b')
% xlabel('Time')
% ylabel('normalization and average I')
% figure(3)
% plot(stept*(1:M),QI,'r',stept*(1:M),PR,'b')


% speed of the wave
speed=(meanI(end)-meanI0)/(Tmax-T0);
text(meanI(end),0.05,sprintf('V=%g',speed))
text(meanI(end),0.04,sprintf('finf=%g',mean(finf)))
% SDs of I and R
meanR=sum(R.*(1:L))/sum(R);
sdI=sqrt(sum(I.*((1:L)-meanI(end)).^2)/sum(I));
sdR=sqrt(sum(R.*((1:L)-meanR).^2)/sum(R)); 
disp([cz,speed,mean(finf)])
% Speed=[Speed,speed]; 
% SDI=[SDI,sdI];
% SDR=[SDR,sdR];
% Finfec=[Finfec,finf(M)];

