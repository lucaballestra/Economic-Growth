function [lbar,kbar,ybar] = equilibrium(alfa,s,delta,n1,n2)

% computes the steady state as in the standard Solow model

lbar = n1/n2;

ka = 0.00000001;
kb = 1000;

fa = s*(lbar^alfa)*(ka^(1-alfa)) - (delta)*ka;
fb = s*(lbar^alfa)*(kb^(1-alfa)) - (delta)*kb;

for i=1:1:1000
    kc = (ka+kb)/2;
    fc = s*(lbar^alfa)*(kc^(1-alfa)) - (delta)*kc;
    if ((fa*fc)<0)
        kb = kc;
        fb = fc;
    else
        ka = kc;
        fa = fc;
    end;
end;
        
    
kbar = kc;
ybar = (lbar^alfa)*(kbar^(1-alfa));



end

