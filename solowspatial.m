clear all;

% This Matlab script solves a spatial Solow model (along the unit circle) 
% with economic-driven migration of labor and capital. The solution of the 
% model is obtained by numerical approximation (finite difference 
% approximation in both space and time).
% 
% Three different test-cases can be considered, depending on the value of
% the dummy variable TESTCASE (to be set by the user): 

% If TESTCASE = 0 the initial distributions of labor and capital are as follows: 
%
% L = Lbar/2      K = Kbar(1 + 0.25 * cos(theta))
%
% where Lbar and Kbar are the steady state solutions for labor and capital in the
% standard (non-spatial) Solow model;
%
% If TESTCASE = 1 the initial distributions of labor and capital are those 
% of an economy with a developed region, where capital and labor are at the 
% equilibrium level, and a capital-rich and scarcely populated region, where 
% capital is at the equilibrium level and labor is significantly below the 
% equilibrium level; 

% If TESTCASE = 2 the initial distributions of labor and capital are those 
% of an economy with both a developed region and an underdeveloped region, 
% where capital and labor are significantly below the equilibrium level.

% the output variables are:
%
% lmat: labor in space and time
% kmat: capitale in space and time
% wagemat: wage in space and time
% kmat: return on investments in space and time
% Lint: total labor in the whole economy (integrated in space) 
% yint: total production in the whole economy (integrated in space) 

TESTCASE  = 0;

% Set the model parameters:
T = 1500; % time horizon
aret = 1; % capital's propensity to migrate 
b = 1; % labor's propensity to migrate
s = 0.3; % saving rate
delta = 0.2; % depreciation rate
alfa = 1/3; % Cobb Douglas coefficient
n1 = 0.02; % coefficient of the logistic population growth dL/dt = n1 * L - n2 * L^2
n2 = 0.01; % coefficient of the logistic population growth dL/dt = n1 * L - n2 * L^2

Nx = 26 +1; % number of discretization space points
Nt = 500000; % number of discretization time levels
Nt0 = Nt;

d = 2*pi;

% computes the space and time discretization steps
dx = d/(Nx-1);
dt = T/(Nt-1);
tvet0 = (T/(Nt0-1))*[0:1:Nt0-1];

for i=1:1:Nx
    xvet(i) = dx*(i-1);
end;
xvet0 = xvet;

% find the (spatially homogeneous) equilibrium of the Solow model
[lbar,kbar,ybar] = equilibrium(alfa,s,delta,n1,n2);

lmat = zeros(Nt,Nx);
kmat = zeros(Nt,Nx);
wagemat = zeros(Nt,Nx);
retmat = zeros(Nt,Nx);

% sets the initial distribution of labor and capital (depending on the
% value of the dummy variable TESTCASE

if (TESTCASE == 0)
k0 = kbar + kbar*cos(xvet)/4;
l0 = 0*k0 + lbar/2;
end;
if (TESTCASE == 1)
    xx1 = pi/4;
ddx = pi/2;
l1 = 1*lbar;
l2 = 0.1*lbar;
k1 = 1*kbar;
k2 = 1*kbar;
k0 = mollifier(xx1,ddx,k1,k2,d,xvet);
l0 = mollifier(xx1,ddx,l1,l2,d,xvet);
end;
if (TESTCASE == 2)
xx1 = pi/4;
ddx = pi/2;
l1 = 1.0*lbar;
l2 = 0.1*lbar;
k1 = kbar;
k2 = kbar/20;
k0 = mollifier(xx1,ddx,k1,k2,d,xvet);
l0 = mollifier(xx1,ddx,l1,l2,d,xvet);
end;

% computes the initial production, wage and return on investments
y0 = (l0.^alfa) .* (k0.^(1-alfa)); 
wage0 = (alfa) * (l0.^(alfa-1)) .* (k0.^(1-alfa)); 
ret0 = (1-alfa) * (l0.^(alfa)) .* (k0.^(-alfa)); 
k = k0;
l = l0;

kold = k;
lold = l;
    
    kt = 1;
    lmat(kt,:) = l0;
    kmat(kt,:) = k0;
    wagemat(kt,:) = wage0;
    retmat(kt,:) = ret0;
    yint(kt) = dx*(2*sum(y0(2:2:Nx-1))/3 + sum(y0(3:2:Nx-2))/3 + (y0(1) + y0(Nx))/6);
    Lint(kt) = dx*(2*sum(l0(2:2:Nx-1))/3 + sum(l0(3:2:Nx-2))/3 + (l0(1) + l0(Nx))/6);
 
% performs the explicit time discretization procedure 
% (finite difference approximation)

for kt = 2:1:Nt
    wage = alfa * (lold.^(alfa-1)) .* (kold.^(1-alfa));
    ret = (1-alfa) * (kold.^(-alfa)) .* (lold.^(alfa));
    yold = (lold.^alfa) .* (kold.^(1-alfa)); 
    % computes the flux of labor 
    for i = 2:1:Nx-1
        derwdx(i) = (wage(i+1)-wage(i-1))/(2*dx); 
    end;
    derwdx(1) = (wage(2)-wage(1))/(dx);
    derwdx(Nx) = (wage(Nx)-wage(Nx-1))/(dx);
    
    % computes the flux of capital
    for i = 2:1:Nx-1
        derfluxdx(i) = (derwdx(i+1)-derwdx(i-1))/(2*dx);
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
        k(i) = (kold(i) + (s*yold(i))*dt -aret*derretdx(i)*dt)/(1+delta*dt);
    end;
    
    l(1) = (l(2) + l(Nx-1))/2;
    k(1) = (k(2) + k(Nx-1))/2;
    l(Nx) = l(1);
    k(Nx) = k(1);
    
    y = (l.^alfa) .* (k.^(1-alfa)); 
    yint(kt) = dx*(2*sum(y(2:2:Nx-1))/3 + sum(y(3:2:Nx-2))/3 + (y(1) + y(Nx))/6);
    Lint(kt) = dx*(2*sum(l(2:2:Nx-1))/3 + sum(l(3:2:Nx-2))/3 + (l(1) + l(Nx))/6);
    
    
    lold = l;
    kold = k;
    
    lmat(kt,:) = l;
    kmat(kt,:) = k;
    wagemat(kt,:) = wage;
    retmat(kt,:) = ret; 
    end;



%lmat = zeros(Nt0,Nx);
%kmat = zeros(Nt0,Nx);
%wagemat = zeros(Nt0,Nx);



%for kt = 1:1:Nt0
%    t = tvet0(kt);
%    lmat(kt,:) = lmat(kt,:);
%    kmat(kt,:) = kmat(kt,:);
%    wagemat(kt,:) = wagemat(kt,:);
%    retmat(kt,:) = retmat(kt,:);
    
%end;

% plots the computed macroeconomic quantities
buildfigures

