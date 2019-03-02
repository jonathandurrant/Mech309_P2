function  [th,r,R_g,V_g] = solve_two_body_prob(t,a,e,p,t0,Cgp);
% This is where you'll sove for R and V (in Frame g) at each time step.

orbital_constants
Re = cst.Re;
mu1 = cst.mu1;

% For now, make something up
th = t/( 2*pi*sqrt(a^3/mu1) )*2*pi;
r = p/(1 + e*cos(th));

% Compute R and V in frame ``p", the perifocal frame.
R_p = [r*cos(th) r*sin(th) 0].'; % m
V_p = [-sin(th) cos(th) + e 0].'*sqrt(mu1/p); % m/s

% Compute R and V in frame ``g", the ECI frame. 
R_g = Cgp*R_p; % m
V_g = Cgp*V_p; % m/s
