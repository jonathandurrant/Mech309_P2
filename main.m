%% MECH 309 Mini-Project 2

close all
clear all
clc

% Format of output
% format short
% format long

%% orbital constants
orbital_constants

%% Load Data
load('MECH309_MP2_data.mat');

whos


%% Solve for velocity given two positions at two times
[rg1,vg1] = find_v_given_position_data(r_g_at_t1,r_g_at_t2,t1,t2);
 

%% Solve for orbital elements
% Note, this function ``orbital_elements" is NOT the same as in MP1; this
% one outputs Delta_t, and not t0.
[a,e,Omega,inc,omega_orbit,Delta_t0] = orbital_elements(rg1,vg1);

% To find t0, need to subtract the time the first radar meas was taken from
% Delta_t0
t0 = Delta_t0 - t1;

% Period
T = 2*pi*sqrt(a^3/cst.mu1);

%% Solve for position given GPS measurments
lv0 = 1; % extra counter.
lv1 = 1;
while (lv0 <= 100)&&(lv1 <= length(t)) % length(t)
    
    % Initial guess from orbit propagation
    [R_orbit,V_orbit] = orbit_propagation(a,e,Omega,inc,omega_orbit,t0,t(lv1));    
    
    % Iterate to find orbit data for each row in each measurement matrix
    
    % Extract measurment data
     meas_data = measurements(:,:,lv1);
     meas_data_size = size(meas_data);
     N_meas = meas_data_size(1,1);
    
    % R_g at time 1
    % r_g_t1_corrected = nonLinearLS(8,0.001);
    % R_g at time 1
    % r_g_t2_corrected = nonLinearLS(9,0.001);
  
   
    
    % Made up, need to change
    SC_r_g_initial_hat(lv1,:) = R_orbit; % Initial estimate of the receiver position
    SC_r_g_hat(lv1,:) =  nonLinearLS(lv1,0.01); % Estimate of the receiver position
    bias_hat(lv1) = 10;
    b_error(lv1) = bias_hat(lv1)/299792458; % Receiver bias divided by speed of light
    
    s
    
    % Update counters
    lv1 = lv1 + 1; 
    lv0 = lv0 + 1;
end

plot_script_v1;

