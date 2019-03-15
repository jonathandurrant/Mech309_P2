% GREG

function root = secantE (M, e)
%% Solving for root E using Secant Method

% Function definition 
f = inline('x - a*sin(x) - b','x', 'a', 'b');

% Initial guesses (in rad): 0 and 90 degrees 
x(1)= 0;
x(2)= 1.57;

%% Start Iterations
tol = 10^(-2);
iteration = 0;
max_iter = 1000;

for lv2=3:max_iter
   x(lv2) = x(lv2-1) - (f(x(lv2-1), e, M))*((x(lv2-1) - x(lv2-2))/(f(x(lv2-1), e, M) - f(x(lv2-2), e, M)));
    iteration=iteration+1;
    abs((x(lv2)-x(lv2-1))/x(lv2))*100
    if abs((x(lv2)-x(lv2-1))/x(lv2))*100<tol
        root = x(lv2);
        iteration = iteration;
        break
    else disp('error');
    end
    
% Return

end

