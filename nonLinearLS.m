function [r_g] = nonLinearLS(matrix_num,tol)
%NONLINEARLS Summary of this function goes here
%   Detailed explanation goes here

%% orbital constants
orbital_constants

%% Load Data
load('MECH309_MP2_data.mat');

whos
%% Build Matrices
for i = 1:6
    
  % extract orbital elements out of the row
    rho(i) = measurements(i,1,matrix_num);
    sigma(i) = measurements(i,2,matrix_num);
    a = measurements(i,3,matrix_num);
    e = measurements(i,4,matrix_num);
    Omega = measurements(i,5,matrix_num);
    inc = measurements(i,6,matrix_num);
    omega_orbit = measurements(i,7,matrix_num);
    t0(i) = measurements(i,8,matrix_num);
    b_i(i) = measurements(i,9,matrix_num);
    
    % Generate position vector for each row
    [R_orbit,V_orbit] = orbit_propagation(a,e,Omega,inc,omega_orbit,t0(i),t(matrix_num));
    
    % extract ith satelite position data (row by row)
    X(i) = R_orbit(1);
    Y(i) = R_orbit(2);
    Z(i) = R_orbit(3);
    
end 

%% Build jacobian and b matrix 

% first guess at spacecraft reciever variables 
x_r = r_g_at_t1(1);
y_r = r_g_at_t1(2);
z_r = r_g_at_t1(3);
% b_t = 0.1;

delta_x = 1; 
iter = 1; 
tol = 0.1;

while iter < 50 && norm(delta_x) > tol
    for i = 1:6

        % calculate distance between ith satelite and reciever 
        d_i = sqrt( (x_r - X(i))^2 + (y_r - Y(i))^2 + (z_r - Z(i))^2);
    
        % Jacobian Matrix - Populate one row at a time
        A_i(i, 1) = (x_r - X(i))/d_i;
        A_i(i, 2) = (y_r - Y(i))/d_i;
        A_i(i, 3) = (z_r - Z(i))/d_i;
        % A_i(i, 4) = 1;

        % calculate b( variable: B) (for non linear least squares: f(x) = b)
        B(i,1) = rho(i) + b_i(i) - sigma(i);
        RHS(i,1) = B(i,1) - d_i; % - b_t;

    end 
    
    % Calculate delta x variable
    delta_x = (A_i'*A_i)\(A_i'*RHS);

    %update variables
    x_r = x_r + delta_x(1);
    y_r = y_r + delta_x(2);
    z_r = z_r + delta_x(3);
    % b_t = b_t + delta_x(4);
    iter = iter+1;
    
end
t = t(matrix_num);
r_g = [A_i(i, 1); A_i(i, 2); A_i(i, 3)]; 
  

end 