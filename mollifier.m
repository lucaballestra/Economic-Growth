function [y] = mollifier(x1,dx,y1,y2,d,xvet)

% builds a C1 function that is equal to:
% y1 in [0, x1] and [d - x1, d]
% y2 in [x1+dx, d-x1-dx]

N = max(size(xvet));
for i=1:1:N
    x = xvet(i);
    if (x<=x1)
        y(i) = y1;
    end;
    if (x>x1) && (x<(x1+dx))
        y(i) = y1 + (y2-y1)*(sin((x-x1)*pi/(2*dx)))^4; 
    end;
    if (x>=(x1+dx)) && (x<=(d-x1-dx))
        y(i) = y2;
    end;
    if (x>(d-x1-dx)) && (x<(d-x1))
        y(i) = y1 - (y1-y2)*(sin((x-(d-x1))*pi/(2*dx)))^4; 
    end;
    if (x>=(d-x1))
        y(i) = y1;
    end;

end

