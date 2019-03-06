function [rg1,vg1] = find_v_given_position_data(rg1,rg2,t1,t2)
% Find vg1 given rg1, rg2, t1, and t2.
% Once rg1 and vg1 are found, the orbital elements can be computed. 

orbital_constants
load('MECH309_MP2_data.mat');
Re = cst.Re;
mu1 = cst.mu1;


% Here is where you will write code to find vg1 given rg1, rg2, t1, and t2.

%% Previous version

% Calculate 2-norm for rg1 and rg2
rg1 = norm(r_g_at_t1);
rg2 = norm(r_g_at_t2);

% Theta in RAD = angle between vectors calculated - via dot product 
% Change in time between the two points 

time_elapse = t2 - t1;
theta = acos((dot(r_g_at_t1, r_g_at_t2)) / (rg1 * rg2));

% magnitude of rg1 X rg2 (cross product)
cross_mag =  rg1 * rg2 * sin(theta);

% Value of T
T = 0.5 * cross_mag;

%
% Two-body problem closed formula
m = (mu1 * (time_elapse)^2) / ( (2*sqrt( rg1 * rg2) * cos(theta / 2))^3 ); 
l = (rg1 + rg2) / ( 4 * (rg1 * rg2)^0.5 * cos(theta/2) ) - 0.5;


% calculate eta using secant method with tolerance = tol
tol = 10 ^-8; 

% eta should be a positive real number. 

eta = secant(1, 10, tol , m, l, 1, 1000);

%test the eta result by using the bisection method. 
% bisectionComparison = bisection(1, 1.0273, 0.00001, m, l);

etaRootTest = eta_f(eta,m,l)/ eta

% semilatus rectum (p) 

p = (eta * cross_mag)^2 / (mu1 * time_elapse^2)

% calculate F and G

F_num = p - rg2 * (1 - cos(theta)); 
F = F_num / p;

G_num = rg1 * rg2 * sin(theta);
G_denom = sqrt (mu1 * p);
G = G_num / G_denom;

% Compute vg1 and output rg1 vector. 

vg1 = (r_g_at_t2 - F * r_g_at_t1) / G
rg1 = r_g_at_t1;

% plot 

% x = 0.5:1.015:1.1;
% plot(x, eta_f(x,m,l));
% grid on
end

