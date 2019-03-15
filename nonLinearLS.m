function [r_g, b_r, iter] = nonLinearLS(matrix_num,tol, rg_init)
%NONLINEARLS Summary of this function goes here
%   Detailed explanation goes here

%% orbital constants
orbital_constants

%% Load Data
load('MECH309_MP2_data.mat');


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
x_r = rg_init(1);
y_r = rg_init(2);
z_r = rg_init(3);
b_r = tol +1;

delta_x = 1; 
iter = 1; 


while iter < 100 && norm(delta_x) > tol
    for i = 1:6

        % calculate distance between ith satelite and reciever 
        d_i = sqrt( (x_r - X(i))^2 + (y_r - Y(i))^2 + (z_r - Z(i))^2);
     
        % Jacobian Matrix - Populate one row at a time
        A_i(i, 1) = (x_r - X(i))/d_i;
        A_i(i, 2) = (y_r - Y(i))/d_i;
        A_i(i, 3) = (z_r - Z(i))/d_i;
        A_i(i, 4) = 1;

        % calculate b( variable: B) (for non linear least squares: f(x) = b)
        B(i,1) = rho(i) + b_i(i) - sigma(i);
        RHS(i,1) = B(i,1) - d_i - b_r;

    end 
    
    % Calculate delta x variable
    delta_x = (A_i'*A_i)\(A_i'*RHS);

    %update variables
    x_r = x_r + delta_x(1);
    y_r = y_r + delta_x(2);
    z_r = z_r + delta_x(3);
    b_r = b_r + delta_x(4);
    iter = iter+1;
    
    
    
end

t = t(matrix_num);
r_g = [x_r, y_r, z_r]; 
  

end 