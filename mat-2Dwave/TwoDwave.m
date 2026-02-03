% Model of escape wave in anisotropic plane xy: x direction of escape,
% Based on 1D model Lin et al JTB 2003
% Mutations in discrete genetic space 
% x fluctuates around lattice sites
%% Model parameters: D=1e-2, a=15, N=1e8 no ordering required in 1D
R0=2.6;                  % basic reproduction ratio
D=0.0001;                % mutation rate
a=7;                  % immunity half-distance
N=1e10;                 % total pop.size for  cutoff of tails
varx=0.01;                 % max amplitude of variation 
%asymmetry='ani.asy';   % anysotropic, asymmetric in x
%asymmetry='ani.sym';   % anysotropic, symmetric in x
asymmetry='iso';       % isotropic: 2 antigenic coordinates
 
%% Internal parameters
L=50;                   % number of variants for each coordinate
Tmax=600;               % max time. [unit] = 1/nu = 1 wk
tshow=100;              % SHOW AT INTERVALS
T0=00;                 % start plotting
stept=1;                % step in time 
M=Tmax/stept;           % number of time points
SEED = 3;
%% Immunity matrix
rng(SEED); % ошибка трансяционной репрессии
%% Projection of lattice on x + disorder
X=(ones(L,1)*(1:L)+varx*(2*rand(L,L)-1));
% Projection on y + disorder
Y=((1:L)'*ones(1,L)+varx*(2*rand(L,L)-1)); 
%% Projection of lattice onto diagonal + disorder
%X=((ones(L,1)*(1:L)+(1:L)'*ones(1,L))/sqrt(2)+varx*(2*rand(L,L)-1));
% Projection onto another diagonal + disorder
%Y=((ones(L,1)*(1:L)-(1:L)'*ones(1,L))/sqrt(2)+varx*(2*rand(L,L)-1));
%%
% Immunity 4D-matrix
K=zeros(L,L,L,L);
for i=1:L
    for j=1:L
        switch asymmetry
            case 'ani.sym'
                dist=abs(X-X(i,j))/a;
                K(i,j,:,:)=dist./(1+dist);
            case 'ani.asy'
                dist=abs(X-X(i,j))/a;
                K(i,j,:,:)=(X < X(i,j)).*dist./(1+dist);
            case 'iso'
                dist=sqrt((X-X(i,j)).^2+(Y-Y(i,j)).^2)/a;
                K(i,j,:,:)=dist./(1+dist);
        end
    end
end

%%
% Initial conditions
% %% near left border
  I=zeros(L,L); S=zeros(L,L); 
%  I(L/2,10)=1e-2;               % Infected dot
%  S(L/2,1:9)=(1-1e-2)/10;      % Susceptible tail
%   I(:,16)=1e-2/L;               % Infected line
%   S(:,1:15)=(1-1e-2)/15/L;      % Susceptible strip
%% elipsis in center
S=zeros(L,L); 
I=((X-L/2).^2+(Y-L/2).^2 < 5^2); I=I/sum(sum(I));
%% near bottom left corner
% I = (X > 25/a & X < 27/a); I = 1e-2*I/sum(sum(I)); % Infected front
% S = X < 25/a; S = S*0.99/sum(sum(S));              % Susceptible corner
%%
R=S;
Inew=zeros(L,L);
norm=zeros(1,M);
finf=zeros(1,M);
% figure(1);pcolor(reshape(K(L/2,L/2,:,:),L,L));colorbar
% figure(2);pcolor(reshape(Kiso(L/2,L/2,:,:),L,L));colorbar
figure(20);
xlabel('Immunity matrix')
subplot(2,2,1)
pcolor(S);colorbar
xlabel('Succeptible')
subplot(2,2,2)
pcolor(I);colorbar
xlabel('Infected')
subplot(2,2,3)
pcolor(R);colorbar
xlabel('Recovered')

%% Main calculation
% Loop in time
for k=1:M
    norm(k)=sum(sum(I+S));          % should be = 1
    finf(k)=sum(sum(I))/norm(k);    % fraction of infected
    % Loop in x and y
    for i=1:L
        for j=1:L 
            Q=sum(sum(S.*reshape(K(i,j,:,:),L,L)));
            P=sum(sum(I.*reshape(K(:,:,i,j),L,L)));
            Inew(i,j)=I(i,j)*(1+stept*(R0*Q-1)); 
            S(i,j)=S(i,j)*(1-stept*R0*P)+stept*I(i,j);
            R(i,j)=R(i,j)+stept*I(i,j);
			% Mutation term
            if j > 1 && j < L && i > 1 && i < L
                Inew(i,j)=Inew(i,j)+stept*D*(I(i,j+1)+I(i,j-1)+I(i-1,j)+I(i+1,j)-4*I(i,j)); %+I(i,j-1)+I(i-1,j)+I(i+1,j)-3*I(i,j));
            end
        end
    end % loop in x and y
    % cutoff of each strain at 1 copy
    Inew=Inew.*(Inew/norm(k) > 1/N);
    % update infected
    I=Inew;  
    % plot in intervals of time tshow after T0
    num=stept*k/tshow;
    if round(num)==num & stept*k > T0;
       figure(20+num)
            subplot(2,2,1)
        pcolor(S)
        colorbar
        xlabel('Susceptible')
        title(sprintf('R0=%g D=%g a=%g time=%g \n norm=%4.2f finf=%6.6f N=%g varx=%g as=%s',R0,D,a,k*stept,norm(k),finf(k),N,varx,asymmetry))
            subplot(2,2,2)
        pcolor(I)
        colorbar
        %title(sprintf('R0=%g D=%g a=%g T=%g %g] norm=%4.2f finf=%6.6f N=%g varx=%g as=%s',R0,D,a,T0,Tmax,norm(i),finf(i),N,varx,asymmetry))
        xlabel('Infected')   
            subplot(2,2,3)
        pcolor(R)
        colorbar
        xlabel('Recovered')
        %title(sprintf('R0=%g D=%g a=%g T=[%g %g] norm=%4.2f finf=%6.4f N=%g varx=%g as=%s',R0,D,a,T0,Tmax,norm(i),finf(i),N,varx,asymmetry))
    end
end  % loop in time


