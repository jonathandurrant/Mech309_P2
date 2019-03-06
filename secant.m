function outputArg1 = secant(x_0,x_1, tol, m, l, Case, iter)
%SECANT Approximates the roots using the secant method
%   x_2 = x_n, x_1 = x_(n-1), x_0 = x_(n-2)

%% From MP2
if ( Case == 1 ) &&  ( eta_f(x_1,m,l)*eta_f(x_0,m,l) < 0)
    x_2 = x_1 - ( eta_f(x_1,m,l) * (x_1 - x_0) )/ ( eta_f(x_1, m,l) - eta_f(x_0,m,l) );
    n = 1;
    
    % check for existence

    while (n < iter) && (abs(x_2 - x_1) > tol)
        x_0 = x_1;
        x_1 = x_2;
        x_2 = x_1 - ( eta_f(x_1,m,l) * (x_1 - x_0) )/ ( eta_f(x_1,m,l) - eta_f(x_0,m,l) );
        n = n+1;
    end

    outputArg1 = x_2;
    iter = n;
end 



end

