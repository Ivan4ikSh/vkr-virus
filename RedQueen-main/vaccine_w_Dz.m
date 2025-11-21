% Red Queen vaccine effect, last vaccination only, a finite width of peak w><bz,
% vaccine delay Dz << bz, w (S6 Txt)

%% Parameters
R0=1.8;             % basic reproduction number without vaccine, H3N2
ez=0.4;             % vaccine max. efficacy
py=0.024;           % susceptivility to homologous strain
by= 1.7;%4.1;%6.4; %1.7;%2.0;  % K(u) scale, best-fit for different g(s)
bz=1.7;
Ub=10^(-4.5);%1e-3;%(-2.6);%(-4.5);       % ben. mutation rate, best fit for different g(s)
N=1e8;
U=Ub;               % total mutation rate in asymmetry \beta
distribution_s='constant';%'uniform';%'constant';%'halfgaussian';

w20=2;            % initial value SD of infected peak, the best-fit without vaccine, 
du=0.1;    
u=0:du:200; % antigenic coordinates (denoted -u)
% natural immunity transmission rate

K=(1+(1-py)/py*exp(-u/by)).^(-1);   % eq. 8
Kder=(1-py)/py/by*exp(-u/by).*(1+(1-py)/py*exp(-u/by)).^(-2);  
cz=0:0.02:1; %  coverage 

%% loop in vaccine breadth b_z
kcol=0;
col='krmbkrmb';
for Dz=[0.5 0.2 0.1 0] 
    kcol=kcol+1;
    Kz=(1+ez/(1-ez)*exp(-u/bz)).^(-1);  % eq. 10
    Kzder=ez/(1-ez)/bz*exp(-u/bz).*(1+ez/(1-ez)*exp(-u/bz)).^(-2);  

%% Iterating in w
    w2=w20*ones(1,length(cz)); % initial
    w2old=10*w2;
    Re=ones(1,length(cz)); 
    iw=0;
    while max(abs(w2old./w2-1)) > 1e-4 && iw < 100
        iw=iw+1;
% Reproduction number
        for k=1:length(cz)
            Re(k)=R0*(1-cz(k)*sqrt(2/pi/w2(k))*du*sum(exp(-u.^2/2/w2(k)).*(1-Kz)));  % eq. S37
        end
%% Normalization of the recovered,
        AA=zeros(size(cz));
        sigma1=zeros(size(cz));
        for k=1:length(cz)
% Recovered individual density  
            Anew=1; A=77; i3=0;
            while abs(A/Anew-1)>1e-4 && i3<200
                i3=i3+1; A=Anew;
                r=exp(-A*Re(k)*du*cumsum(K));  % eq. S5, denoted s
                totalr=sum(r)*du; 
                Anew=1/totalr;
            end
            r=A*r;
            AA(k)=A; 
% Selection coefficient from natural infections but with reduced Re
            sigma1(k)=Re(k)*du*sum(Kder.*r); % eq. S7, first term
        end % loop in cz
        sigma0=sigma1(1);   % sigma without vaccine, cz=0
        A0=AA(1);           % normalization factor without vaccine
        
%% Finding full sigma including vaccination
        xi=sqrt(2/pi)*Dz*(w2).^(-0.5);
        eta=sqrt(2/pi)*bz/(1-ez)*(w2).^(-0.5);
        Re_der=R0*ez*(1-ez)/bz*cz.*min([ones(1,length(cz)); xi; xi.*eta]); % asymptotics Eqs S22, S44, S46 stitched together
        sigma=sigma1 + Re_der./Re; % eq. S7 

%% Update w2
        w2old=w2;
% Initial interation of speed and fraction of infected finf
        V=2*sigma*log(N); 
        finf=V.*AA;
% Iterating finf
        fold=0.1*finf; 
        i_f=0;
        while max(abs(finf/fold-1)) > 1e-3 && i_f<100
            i_f=i_f+1;
            fold=finf;
            Ninf=N*finf;
            switch distribution_s
                case 'constant'
                    for i1=1:20 % iterating in V
                        V=2*sigma.*log(Ninf./sqrt(V.^2 .*log(V/Ub)/(sigma.^3*Ub)))./(log(V/exp(1)/Ub).^2+1); % if V > sigma
                    end
                    w2=V./sigma; % variance in antigenic coordinate
                    s_ast=sigma;
                case 'exponential' 
                    v=2*sigma.^2 .*log(Ninf*Ub);     % Initial iteration of adaptation speed
                    s_ast = sqrt(2*v.*log(2/Ub*sqrt(v/2/pi)));  % Most probable fixed allele selection coefficient
                    for i1=1:10
                        xc = s_ast+v./sigma;
                        s_ast = sqrt(2*v.*log(2/Ub*sqrt(v/2/pi)/(1+sigma./xc+v./sigma./s_ast)));
                        v=2*sigma.*(-s_ast+sigma.*log(N*Ub*xc.^2 ./v-1+2*xc.*sigma./v+2*sigma.^2 ./v));
                    end
                    V=v./s_ast;
                    w2=v./s_ast.^2;    % variance in antigenic coordinate
                case 'halfgaussian'
                    s_ast=sigma*sqrt(pi).*log(sigma/Ub).^0.5;   % Initial iteration
                    kappa=log(2*Ninf.*s_ast*Ub/pi./sigma);
                    for i1=1:10
                        kappa=log(2*Ninf.*s_ast*Ub.*(1+kappa).^1.5/pi./sigma./kappa)./...
                        log(s_ast/Ub.*sqrt(kappa./(1+kappa))); 
                % rederived on 17.12.2024 for half-Gaussian rho(s) with average s=1 above
                        s_ast=sigma.*sqrt(pi*kappa./(1+kappa)).*log(s_ast/Ub.*sqrt(kappa./(1+kappa))).^0.5;  % rederived on 17.12.2024
                    end
                    v=pi*sigma.^2/2 .*kappa;
                    V=v./s_ast;
                    w2=v./s_ast.^2;    % variance in antigenic coordinate
                case 'uniform'
                    alpha=10; % step steepness
                    s0=sigma*2*gamma(1+1/alpha)/gamma(1+2/alpha); % to have mean(s)=sigma
                    X=sqrt(alpha)*s0/Ub;
                    % adaptation rate
                    v=2*s0.^2 .*log(Ninf.*s0.*sqrt(log(Ninf.*s0)))./log(X.*log(X).^0.5).^2 ./(1+log(X)/4 ./log(Ninf.*s0)).^2;
                    % effective selection coefficient
                    s_ast=s0.*(sqrt(2)/alpha*log(s0/Ub)).^(1/(alpha-1));
                    V=v./s_ast;
                    w2=v./s_ast.^2;    % variance in antigenic coordinate
            end  % switch distribution s
            finf=V.*AA; % updating finf
        end % iterations in finf
    end % iterations in w2
    
%% Plotting results
    ii=1:length(Re);%find(R0 > 1);
figure(1)
    plot(cz(ii), V(ii)/V(1), col(kcol))
    text(cz(round(end)),V(round(end))/V(1),sprintf('%g', Dz))
    hold on
    xlabel('Coverage c_z')
    ylabel('Normalized evolution rate and normalized incidence')
    title(sprintf('R0=%g, e_z=%g, p_y=%g, b_y=%g, b_z=%g, N=%g, U_b=%g \n %s, Dz on curves' ,R0,ez,py,by,bz,N,Ub,distribution_s ))
    axis([0 1 0 3])
    hold on
    plot(cz(ii), finf(ii)/finf(1), [col(kcol) '--'])
    text(cz(round(end)),finf(round(end))/finf(1),sprintf('%g', Dz))
end % loop in Dz
figure(1)
box off
hold off



    

