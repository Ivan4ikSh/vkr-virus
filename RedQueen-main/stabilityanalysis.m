% eigenvalues of red queen problem
Re=1.1;        % effective reproduction number (to be changed in a loop later)
py=0.024;       % susceptibility to homologous strain
by=1.7;%2.0;         % steepness of K(u) 
du=0.02;
u=0:du:50;      % antigenic coordinates 
K=(1+(1-py)/py*exp(-u/by)).^(-1);

% Recovered individual density 
Anew=1;A=77;i=0;
while abs(A/Anew-1)>1e-6 && i<100
    i=i+1;
    A=Anew;
    r=exp(-A*Re*du*cumsum(K));
    totalr=sum(r)*du; 
    Anew=1/totalr;
end
r=A*r; % final recovered density = exp(Phi0)
% derivative log
derPhi0=A*Re*K; 
% Eigenvalues
xmax = 1; % max size of real part of eigenvalue
xmin = -0.6;  % min size of real part
dx = 0.005; 
ymax = 3; % max size of imag part
dy = 0.2; 
x = xmin:dx:xmax;
y= 0:dy:ymax;
LHS=zeros(length(y), length(x));

for k=1:length(y)
    for i=1:length(x)
        lambda = x(i)+sqrt(-1)*y(k);
        int13 = - du*cumsum(exp(lambda*u).*derPhi0);
        LHS(k,i) = du*sum(r.*exp(-lambda*u).*(1+int13));
    end
end

figure(3)
% xRe LHS and Im LHS as a function of lambda
[C,h]=contour(x,y,real(LHS),[1e10 0],'r');
%clabel(C,h)
hold on
[C,h]=contour(x,y,imag(LHS),[1e10 0],'b');
%clabel(C,h)
xlabel('Re(lambda)')
ylabel('Im(lambda)')
title(sprintf('real(LHS) red imag(LHS) blue Re=%g py=%g by=%g',Re,py,by))
grid 
hold off
% 
% figure(4)
% [C,h]=contour(x,y,abs(LHS).^2,[0 1e-5 1e-4 ],'b');
% clabel(C,h)
% xlabel('Re(lambda)')
% ylabel('Im(lambda)')
% title(sprintf('abs(LHS)^2 Re=%g py=%g by=%g',Re,py,by))
% grid 

