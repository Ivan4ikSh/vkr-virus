% vaccine effect, the last vaccination only

%% Parameters
R00=1.8;        % basic reproduction number without vaccine
e_z=0.4;        % vaccine max. efficacy
py=0.024;       % susceptivility to homologous strain
by=1.7;%4.1;%6.8;%1.7;%2.0;         % steepness of K(u) 
vaccineahead=1;
du=0.05;
u=0:du:1000;%10000;      % antigenic coordinates 
% natural immunity transmission rate
%  K=1-exp(-u);  
%  Kder=exp(-u);
K=(1+(1-py)/py*exp(-u/by)).^(-1);
Kder=(1-py)/py/by*exp(-u/by).*(1+(1-py)/py*exp(-u/by)).^(-2);
c_z=0:0.02:1; %  coverage 
R0=R00*(1-e_z*c_z);

%% normalization and 1st component of sigma
AA=zeros(size(c_z));sigma1=zeros(size(c_z));
for k=1:length(c_z)
% Recovered individual density 
    Anew=1;A=77;i=0;
    while abs(A/Anew-1)>1e-6 && i<1000
        i=i+1;
        A=Anew;
        r=exp(-A*R0(k)*du*cumsum(K));
        totalr=sum(r)*du; 
        Anew=1/totalr;
    end
    r=A*r;
    AA(k)=A; 
% Selection coefficient from derivative in K
    sigma1(k)=R0(k)*du*sum(Kder.*r);
    Norm(k)=sum(r)*du;
end
sigma0=sigma1(1);
A0=AA(1);


%% loop in kappa = (1-e_z)/b_z
kcol=0;
col='krgbkrgb';
for kappa=[0 0.1 (1-e_z)/1.7 0.6];      % vaccine cross-reactivity reduction parameter
    kcol=kcol+1;

%% Finding dependence of x = sigma/sigma0 on beta

% calculating X=sigma/sigma0
  
if vaccineahead % if vaccine is ahead of the peak
    X = (sigma1 - e_z*c_z./(1-e_z*c_z)*kappa)/sigma0;
else
    X = (sigma1 + e_z*c_z./(1-e_z*c_z)*kappa)/sigma0;
end
    % relative incidence
    incid = X.*AA/A0; 

    figure(1)
    ii=1:length(R0);%find(R0 > 1);
    plot(c_z(ii), X(ii), [col(kcol),'-'], c_z(ii), incid(ii), [col(kcol),'--'])
    hold on
    text(c_z(round(end/4)),incid(round(end/4)),sprintf('%g',(1-e_z)/kappa))
   ylabel('Normalized substitution rate X (-) Incidence X*A/A_0 (.-) R_e (:) ')
    title(sprintf('K_s Myrthe, R0=%g, e_z=%g, p_y=%g, b_y=%g, b_z on curves',R00,e_z,py,by ))
end % loop in b_z
box off
hold off
axis([0 min((1-1/R00)/e_z,1) 0 3])

xlabel('Coverage c_z')

figure(2)

plot(c_z(ii), R0(ii))
xlabel('Coverage c_z')
ylabel('R_e')
box off



    

