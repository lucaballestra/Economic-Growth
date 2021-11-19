% This script, once the model has been solved numerically, plots
% the space-time distributions of labor, capital, wages and return on
% investments (Figures 1-4) and the time evolution of the total labor
% and total production (Figures 5, 6)


set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

figure(1)
mesh(xvet,tvet0,lmat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$L$');

figure(2)
mesh(xvet,tvet0,kmat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$K$');

figure(3)
mesh(xvet,tvet0,wagemat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$w$');

figure(4)
mesh(xvet,tvet0,retmat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$r$');

figure(5)
plot(tvet0, Lint);
xlabel('$t$');
ylabel('$L_{tot}$');

figure(6)
plot(tvet0, yint);
xlabel('$t$');
ylabel('$Y_{tot}$');
