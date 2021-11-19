clear all;

% this Matlab script solves the spatial Solow model by numerical
% approximation (finite difference approximation in both space and time)

T = 800; % time horizon
aret = 1; % capital's propensity to migrate 
b = 1; % labor's propensity to migrate
s = 0.3; % saving rate
delta = 0.2; % depreciation rate
alfa = 1/3; % Cobb Douglas coefficient
n1 = 0.02; % coefficient of the logistic population growth dL/dt = n1 * L - n2 * L^2
n2 = 0.01; % coefficient of the logistic population growth dL/dt = n1 * L - n2 * L^2

Nx = 26 +1; % number of discretization space points
Nt = 50000; % number of discretization time levels
Nt0 = Nt;


d = 2*pi;
dx = d/(Nx-1);
dt = T/(Nt-1);

for i=1:1:Nx
    xvet(i) = dx*(i-1);
end;
xvet0 = xvet;

[lbar,kbar,ybar] = equilibrium(alfa,s,delta,n1,n2);


lmat1 = zeros(Nt,Nx);
kmat1 = zeros(Nt,Nx);
wagemat1 = zeros(Nt,Nx);
ymat1 = zeros(Nt,Nx);
retmat1 = zeros(Nt,Nx);

% sets the initial distribution of labor and capital
k0 = kbar + kbar*cos(xvet)/4;
l0 = lbar/2 + 0*lbar*cos(xvet);

wage0 = (alfa) * (l0.^(alfa-1)) .* (k0.^(1-alfa)); 
ret0 = (1-alfa) * (l0.^(alfa)) .* (k0.^(-alfa)); 
k = k0;
l = l0;

kold = k;
lold = l;
    
    y = (lold.^alfa) .* (kold.^(1-alfa)); 
    wage = alfa * (lold.^(alfa-1)) .* (kold.^(1-alfa)); 
    ret = (1-alfa) * (lold.^(alfa)) .* (kold.^(-alfa)); 

    kt = 1;
    lmat1(kt,:) = l;
    kmat1(kt,:) = k;
    wagemat1(kt,:) = wage;
    retmat1(kt,:) = ret;
    ymat1(kt,:) = y;
 
% explicit time discretization step (finite difference approximation)

for kt = 2:1:Nt
    y = (lold.^alfa) .* (kold.^(1-alfa)); 
    wage = alfa * (lold.^(alfa-1)) .* (kold.^(1-alfa));
    ret = (1-alfa) * (kold.^(-alfa)) .* (lold.^(alfa));
    for i = 2:1:Nx-1
        derwdx(i) = (wage(i+1)-wage(i-1))/(2*dx); 
    end;
    derwdx(1) = (wage(2)-wage(1))/(dx);
    derwdx(Nx) = (wage(Nx)-wage(Nx-1))/(dx);
    
    for i = 2:1:Nx-1
        derfluxdx(i) = (derwdx(i+1)-derwdx(i-1))/(2*dx);
        %der2kdx(i) = (kold(i+1)-2*kold(i)+kold(i-1))/(dx*dx);
        
        lp = (lold(i+1)+lold(i))/2;
        lm = (lold(i-1)+lold(i))/2;
        kp = (kold(i+1)+kold(i))/2;
        km = (kold(i-1)+kold(i))/2;
        derfluxdx(i) = ((wage(i+1)-wage(i)) - (wage(i)-wage(i-1)))/(dx*dx);
        derretdx(i) = ((ret(i+1)-ret(i)) - (ret(i)-ret(i-1)))/(dx*dx);
        
        
        
        kp = (kold(i+1)+kold(i))/2;
        km = (kold(i-1)+kold(i))/2;
        der2kdx(i) = (kold(i+1)-2*kold(i)+kold(i-1))/(dx*dx);
        derkdxp = (kold(i+1)-kold(i))/dx;
        derkdxm = (kold(i)-kold(i-1))/dx;
        der2kdx(i) = (derkdxp  - derkdxm )/dx;
        derkwdx(i) = ((wage(i+1)-wage(i)) - (wage(i)-wage(i-1)))/(dx*dx);
    end;
    for  i= 2:1:Nx-1
        x = xvet(i);
        l(i) = lold(i) - b*derfluxdx(i)*dt + dt*(n1*lold(i)-n2*lold(i)*lold(i));
        k(i) = (kold(i) + (s*y(i))*dt -aret*derretdx(i)*dt)/(1+delta*dt);
    end;
    l(1) = (l(2) + l(Nx-1))/2;
    k(1) = (k(2) + k(Nx-1))/2;;
    l(Nx) = l(1);
    k(Nx) = k(1);
    
    lold = l;
    kold = k;
    
    lmat1(kt,:) = l;
    kmat1(kt,:) = k;
    wagemat1(kt,:) = wage;
    retmat1(kt,:) = ret;
    ymat1(kt,:) = y;
    
end;



lmat = zeros(Nt0,Nx);
kmat = zeros(Nt0,Nx);
wagemat = zeros(Nt0,Nx);
ymat = zeros(Nt0,Nx);
lpiccolomat = zeros(Nt0,Nx);
kpiccolomat = zeros(Nt0,Nx);
ksul = zeros(Nt0,Nx);
ksuy = zeros(Nt0,Nx);
consumiperworker = zeros(Nt0,Nx);



tvet0 = (T/(Nt0-1))*[0:1:Nt0-1];

for kt = 1:1:Nt0
    t = tvet0(kt);
    lmat(kt,:) = lmat1(kt,:);
    kmat(kt,:) = kmat1(kt,:);
    wagemat(kt,:) = wagemat1(kt,:);
    retmat(kt,:) = retmat1(kt,:);
    ymat(kt,:) = ymat1(kt,:);
    ksul(kt,:) = kmat(kt,:)./lmat(kt,:);
    ksuy(kt,:) = kmat(kt,:)./ymat(kt,:);
    
    lpiccolomat(kt,:) = lmat(kt,:)-lbar;
    kpiccolomat(kt,:) = kmat(kt,:)-kbar;
    
    % computes savings and consumptions
    for j=1:1:max(size(xvet))
        consumiperworker(kt,j) = (1-s)*ymat(kt,j)./lmat(kt,j);
    end;
    
end;


kperworkermat = kmat./lmat;
wageperworkermat = wagemat./lmat;
yperworkermat = ymat./lmat;

% plots the computed macroeconomic quantities

creafigure

