function  [th,r,R_g,V_g] = solve_two_body_prob_hunter(t,a,e,p,t0,Cgp);
% This is where you'll sove for R and V (in Frame g) at each time step.

orbital_constants
Re = cst.Re;
mu1 = cst.mu1;

%% defining all relevant functions as functions of E
sinTh = @(E) sin(E)*sqrt(1-e^2)/(1-e*cos(E));
cosTh = @(E) (cos(E)-e)/(1-e*cos(E));
theta = @(E) atan2(sinTh(E),cosTh(E));
f = @(E)(E-e*sin(E))-(t-t0)*sqrt(mu1/a^3);

%Using the secant root finding method to find true anomaly
%disp('Secant method to evaluate true anomaly')
%disp('Checking root existence')
hi = pi; 
lo = -pi;
tol = 1e-6;

%if root_exists(f, lo, hi)
%disp('Secant method to evaluate eta')
[E,sectantiter] = secant_f(f, lo, hi, tol);
%else
%end

% plug numerical E value into atan2 for correct side.
th = theta(E);

%% Compute R and V in frame ``p", the perifocal frame.
r = p/(1 + e*cos(th));
R_p = [r*cos(th) r*sin(th) 0].'; % m
V_p = [-sin(th) cos(th) + e 0].'*sqrt(mu1/p); % m/s

% Compute R and V in frame ``g", the ECI frame. 
R_g = Cgp*R_p; % m
V_g = Cgp*V_p; % m/s
