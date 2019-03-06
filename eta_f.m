function y = eta_f(x, m, l)
%ETA_F calculates value of the function f(x) = 1- eta + (m/eta^2)*
%W((m/eta^2) - l)

%% From MP1 

%   x is a placeholder for eta.

w = m ./(x.^2) - l; 

    function y = g
        y = 2 .* asin(sqrt(w));
    end

    function y = W
        y = (2.*g - sin(2 .* g)) ./ (sin(g).^3);
    end

y = 1 - x + (m ./ x.^2) .* W;

end

