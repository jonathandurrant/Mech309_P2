%% From Hunter

function [val, iter] = secant_f (f,x1,x2,tol)
maxiter = 50;
x3 = (x1*f(x2) - x2*f(x1))/(f(x2) - f(x1));
iter = 1; 
while abs(f(x3)) > tol
    x1 = x2;
    x2 = x3;
    x3 = (x1*f(x2) - x2*f(x1))/(f(x2) - f(x1));
    
    iter = iter + 1;
    if(iter == maxiter)
        display('WARNING: MAXITER REACHED')
        break;
    end
end
val = x3;

