% Red Queen vaccine effect, last vaccination only, vaccine coordinate right at peak
% for a given width of peak

%% Parameters
R0=1.8;        % basic reproduction number without vaccine
ez=0.4;        % vaccine max. efficacy
py=0.024;       % susceptivility to homologous strain
by=6.8; %1.7;%2.0;  % K(u) scale, best fit for half-Gaussian distribution s
w0=0.72;            % SD of infected peak, approximate from the best-fit without vaccine
Ub=10^(-2.6);       % total mutation rate, best fit for half-Gaussian distribution s                    

du=0.05;
u=0:du:100; %10000;      % antigenic coordinates (denoted -u)
% natural immunity transmission rate

K=(1+(1-py)/py*exp(-u/by)).^(-1);   % eq. 8
Kder=(1-py)/py/by*exp(-u/by).*(1+(1-py)/py*exp(-u/by)).^(-2);  
cz=0:0.02:1; %  coverage 

%% loop in vaccine breadth b_z
kcol=0;
col='krgbkrgb';
%for bz=[200 6 1.7 1] 
%for U=Ub*[1 2 5 10]
bz=1.7;
U=Ub;
for w=w0*[0.8 1 2 5]
    kcol=kcol+1;
    Kz=(1+ez/(1-ez)*exp(-u/bz)).^(-1);  % eq. 10
    Kzder=ez/(1-ez)/bz*exp(-u/bz).*(1+ez/(1-ez)*exp(-u/bz)).^(-2);
% reproduction number
    Re=R0*(1-cz*sqrt(2/pi)/w*du*sum(exp(-u.^2/2/w^2).*(1-Kz)));  % eq. S37

%% recovered normalization, sigma component from immunity to natural infection
    AA=zeros(size(cz));sigma1=zeros(size(cz));
    for k=1:length(cz)
% Recovered individual density  
        Anew=1; A=77; i=0;
        while abs(A/Anew-1)>1e-5 && i<100
            i=i+1; A=Anew;
            r=exp(-A*Re(k)*du*cumsum(K));  % eq. S5, denoted s
            totalr=sum(r)*du; 
            Anew=1/totalr;
        end
        r=A*r;
        AA(k)=A; 
% Selection coefficient from natural infections but with reduced Re
        sigma1(k)=Re(k)*du*sum(Kder.*r); % eq. S7, first term
%     Norm(k)=sum(r)*du;
    end % loop in cz
    sigma0=sigma1(1);   % sigma without vaccine
    A0=AA(1);           % normalization factor without vaccine

%% Finding X = sigma/sigma0 by iteration 
    XX=zeros(1,length(cz));
    beta0=1/w^6*U/6/sigma0; % eq. S41
    for k=1:length(cz)
        Xnew=1; X=77; i=0; 
        while abs(X/Xnew-1) > 1e-5 && i < 100
            X=Xnew; i=i+1;
            beta=beta0/X;
            Re_der=R0*cz(k)*sqrt(2/pi)*beta/w^2*du*sum(exp(-u.^2/2/w^2).*Kzder); % eq. S42
            Xnew = (sigma1(k) + Re_der/Re(k))/sigma0; % eq. S7 
        end % iteration
        XX(k)=Xnew;
    end % loop in cz
% relative incidence
    incid = XX.*AA/A0; 

    figure(1)
    ii=1:length(Re);%find(R0 > 1);
    plot(cz(ii), XX(ii), [col(kcol),'-'], cz(ii), incid(ii), [col(kcol),'--'])
    text(cz(round(3*end/4)),XX(round(3*end/4)),sprintf('%g', w));%bz))
    text(cz(round(3*end/4)),Re(round(3*end/4)),sprintf('%g', w));%bz))
    hold on
    plot(cz(ii), Re(ii), [col(kcol),':'])
    axis([0 1 0 3])
    xlabel('Coverage c_z')
    ylabel('Normalized substitution rate X (-) Incidence X*A/A_0 (.-) R_e (:) ')
%     title(sprintf('R0=%g, e_z=%g, p_y=%g, b_y=%g, w=%g, Ub+Ud=%g \n b_z on curves',R0,ez,py,by,w,U ))
%     title(sprintf('R0=%g, e_z=%g, p_y=%g, b_y=%g, w=%g, bz=%g \n U/Ub on curves',R0,ez,py,by,w,bz ))
     title(sprintf('R0=%g, e_z=%g, p_y=%g, b_y=%g, Ub=%g, bz=%g \n w on curves',R0,ez,py,by,U,bz ))
end % loop 
box off



    

