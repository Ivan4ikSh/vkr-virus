% vaccine effect, frequent vaccination limit with ovelapping effects 

%% Parameters
R00=1.8;        % basic reproduction number without vaccine
du=0.2;
% K=1-exp(-u);    % natural immunity transmission rate, normalized in x and y
% Kder=exp(-u);
% realistic K(u) fit to influenza data by Myrthe
py=0.024;       % susceptivility to homologous strain
by=1.7;         % steepness of K(u) 
alpha=0:0.05:5; % composite vaccine parameter including coverage, efficiency, frequency 
eps=0.005;        % accuracy of x

u=0:du:300;      % antigenic coordinates
K=(1+(1-py)/py*exp(-u/by)).^(-1);
Kder=(1-py)/py/by*exp(-u/by).*(1+(1-py)/py*exp(-u/by)).^(-2);

u1=0:du:2000;      % antigenic coordinates, long
K1=(1+(1-py)/py*exp(-u1/by)).^(-1);
Kder1=(1-py)/py/by*exp(-u1/by).*(1+(1-py)/py*exp(-u1/by)).^(-2);


%% Case alpha=0
% Iterating A
 Anew=1;A=77;i=0;
 while abs(A/Anew-1)>1e-6 && i<100
      i=i+1;
      A=Anew;
      r=exp(-A*R00*du*cumsum(K));
      totalr=sum(r)*du; 
      Anew=1/totalr;
 end
% Recovered individual density 
r=A*r;
% Normalization factor 
A0=A;
% Selection coefficient 
sigma0=R00*du*sum(Kder.*r);

%% loop in gamma
kcol=0;
%col='krbgkrbg';
for gamma=[0 0.5 1 2 4 8];      % vaccine cross-reactivity reduction parameter
%kcol=kcol+1;

%% Finding dependence of x = sigma/sigma0 on beta
X=zeros(size(alpha)); AA=X; error=X; 
RR=2*ones(size(alpha));
for j=1:length(alpha)

% Iterating x
    k=0; x=1; xnew=2;
    while abs(xnew-x) > eps && k < 100
% Decreased repoduction number
        x=xnew; 
        k=k+1;
        Re=R00*exp(-gamma*alpha(j)/x); 
% Recovered individual density
        if j == 1 || (j > 1 && RR(j-1) > 1.4)
            % normalizing r
            Anew=1;A=77;i=0;
            while abs(A/Anew-1)>1e-6 && i<100
                i=i+1;
                A=Anew;
                r=exp(-A*Re*du*cumsum(K));
                totalr=sum(r)*du; 
                Anew=1/totalr;
            end
            r=A*r;
% Iterating x
            xnew = (Re*du*sum(Kder.*r) + alpha(j)/x)/sigma0;
        else 
             % normalizing r
            Anew=1;A=77;i=0;
            while abs(A/Anew-1)>1e-6 && i<1000
                i=i+1;
                A=Anew;
                r=exp(-A*Re*du*cumsum(K1));
                totalr=sum(r)*du; 
                Anew=1/totalr;
            end
            r=A*r;
% Iterating x
            xnew = (Re*du*sum(Kder1.*r) + alpha(j)/x)/sigma0;
        end % if short or long
    end % iterations end
    error(j)=abs(xnew-x);
    X(j) = xnew;
    AA(j) = A;
    RR(j)=Re;
end % loop in alpha


figure(1)
incid = X.*AA/A0; % relative incidence
ii=find(RR > 1);
% plot(alpha(ii), X(ii), [col(kcol),'-'],...
%     alpha(ii), incid(ii), [col(kcol),'--'])
plot(alpha(ii), X(ii), alpha(ii), incid(ii))
text(alpha(ii(end)),incid(ii(end)),sprintf('%g',gamma))
hold on
xlabel('Vaccine parameter, \alpha')
ylabel('Change in substitution rate (-) and incidence (--)')
title(sprintf(' R0=%g  b_y=%g max(u)=%g max(u1)=%g',R00, by, max(u),max(u1)))

figure(2)
plot(alpha(ii),RR(ii))
text(alpha(ii(end)),RR(ii(end)),sprintf('%g',gamma))
hold on
xlabel('\alpha')
ylabel('R_e')
title(sprintf('  R0=%g b_y=%g max(u)=%g max(u1)=%g',R00, by, max(u),max(u1)))
end % loop in gamma

hold off
box off
axis([0 max(alpha) 1 ceil(R00)])
figure(1)
hold off
box off

