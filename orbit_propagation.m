function [R_g,V_g] = orbit_propagation(a,e,Omega,inc,omega_orbit,t0,t)
% Plot the orbit

orbital_constants

Re = cst.Re;
mu1 = cst.mu1;

%% Compute various useful quantities
rp = a*(1 - e); % m, perigee radius
p = rp*(1 + e); % m, semilatus rectum
h = sqrt(p*mu1);

%% Find the DCM from the perifocal frame to the inertial frame
Cpg = eye(3); % Made up, needs to be changed. 
Cgp = Cpg'; % all angels in radians (rad)

% Compute th, r, R, and V
[th,r,R_g,V_g] = solve_two_body_prob(t,a,e,p,t0,Cgp);
    


end

