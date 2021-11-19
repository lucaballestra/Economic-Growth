% This script, once the model has been solved numerically, plots
% the space-time distributions of labor, capital, wages, return on
% investments, capital per capita and consumption per capita


set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');



figure(1)

mesh(xvet,tvet0,lpiccolomat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$\hat{L}$');

figure(2)
mesh(xvet,tvet0,kpiccolomat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$\hat{K}$');

figure(3)
mesh(xvet,tvet0,kmat./lmat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$\frac{K}{L}$');


figure(4)
mesh(xvet,tvet0,consumiperworker);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$\frac{C}{L}$');



figure(5)
mesh(xvet,tvet0,lmat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$L$');



figure(6)
mesh(xvet,tvet0,kmat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$K$');

figure(7)
mesh(xvet,tvet0,wagemat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$w$');

figure(8)
mesh(xvet,tvet0,retmat);
xlim([0 d]);
xlabel('$\theta$');
ylim([0 T]);
ylabel('$t$');
zlabel('$r$');
