
% Fitting for Red Queen ams  2024, no vaccine, s=const model
% Calculates the mean square diff between V and TMRCA and the experimental values
% Uses the model with a fixed s

function SD2 = call_fitting(X)

global distribution_s R0 N TMRCA_exp V_exp du py TMRCA z V finf K u r Kder sigma x0 s_ast

Ub=10^X(1);
by=X(2);

%% Transmission function K(u) from Myrthe's data inference

u=0:du:150;     % antigenic coordinates
beta=10;    % steepness parameter for the case of a uniform s-distribution

K=(1+(1-py)/py*exp(-u/by)).^(-1);
Kder=(1-py)/py/by*exp(-u/by).*(1+(1-py)/py*exp(-u/by)).^(-2);

%% Iterating normalization factor A
Anew=1; A=77; i=0;
while abs(A/Anew-1)>1e-3 && i<100
        i=i+1;
        A=Anew;
        r=exp(-A*R0*du*cumsum(K));
        totalr=sum(r)*du; 
        Anew=1/totalr;
end

% Recovered individual density 
r=A*r;
% Selection coefficient 
sigma=R0*du*sum(Kder.*r);

%% Substitution rate
V=2*sigma*log(N); % initial approximation

% Iterating fraction of infected
fnew=0.1; finf=1;
while abs(finf/fnew-1) > 1e-2
    finf=fnew;
    switch distribution_s
        case 'constant'
            for i=1:3
                V=2*sigma/(log(V/exp(1)/Ub)^2+1)*...
                   log(N*finf*sqrt(sigma^3*Ub/V^2/log(V/Ub))); % if V > sigma
            end
            x0=V/sigma*log(V/exp(1)/Ub); 
            TMRCA=z*sqrt(2*log(N*finf*sigma)/(V*sigma));
            s_ast=sigma;
        case 'exponential' 
            Ninf=N*finf;
% Initial iteration
            v=2*sigma^2*log(Ninf*Ub);     
            s_ast= sqrt(2*v*log(2/Ub*sqrt(v/2/pi)));
            for i=1:10
                xc = s_ast+v/sigma;
                s_ast = sqrt(2*v*log(2/Ub*sqrt(v/2/pi)/(1+sigma/xc+v/sigma/s_ast)));
                v=2*sigma*(-s_ast+sigma*log(N*Ub*xc^2/v-1+2*xc*sigma/v+2*sigma^2/v));
            end
            x0=xc/s_ast;
            V=v/s_ast; % substitution rate
            TMRCA=z*sqrt(2*log(N*finf*sigma)/v);
        case 'halfgaussian'
            % rho(s)=(1/pi/sigma)*exp(-s^2/pi/sigma^2); average s = sigma
            Ninf=N*finf;
 % Initial iteration
            s_ast=sigma*sqrt(pi)*log(sigma/Ub)^0.5;
            kappa=log(2*Ninf*s_ast*Ub/pi/sigma);
            for i=1:10
                kappa=log(2*Ninf*s_ast*Ub*(1+kappa)^1.5/pi/sigma/kappa)/...
                    log(s_ast/Ub*sqrt(kappa/(1+kappa))); 
                % rederived on 17.12.2024 for half-Gaussian rho(s) with average s=1 above
                s_ast=sigma*sqrt(pi*kappa/(1+kappa))*log(s_ast/Ub*sqrt(kappa/(1+kappa)))^0.5;  % rederived on 17.12.2024
            end
            x0=1+kappa;
            v=pi*sigma^2/2*kappa;
            V=v/s_ast;
            TMRCA=z*sqrt(2*log(N*finf*sigma)/v);
        case 'uniform'
            s0=sigma*2*gamma(1+1/beta)/gamma(1+2/beta); % s0/sigma to have mean(s)=sigma
            X=sqrt(beta)*s0/Ub;
            Ninf=N*finf;
            v=2*s0^2*log(Ninf*s0*sqrt(log(Ninf*s0)))/log(X*log(X)^0.5)^2/(1+log(X)/4/log(Ninf*s0))^2;
            s_ast=s0*(sqrt(2)/beta*log(s0/Ub))^(1/(beta-1));
            V=v/s_ast;
            TMRCA=z*sqrt(2*log(N*finf*sigma)/v);
            x0=77; % no x0 in this case
    end
% New infected fraction
    fnew=V*A;
end

SD2 = sqrt(((TMRCA/TMRCA_exp -1)^2 + (V/V_exp -1)^2)/2);

% if imag(V)>0
%     disp('fuck V') 
% elseif imag(TMRCA)>0
%     disp('fuck TMRCA')
% end



