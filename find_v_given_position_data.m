function [rg1,vg1] = find_v_given_position_data(rg1,rg2,t1,t2)
% Find vg1 given rg1, rg2, t1, and t2.
% Once rg1 and vg1 are found, the orbital elements can be computed. 

orbital_constants

Re = cst.Re;
mu1 = cst.mu1;

% Here is where you will write code to find vg1 given rg1, rg2, t1, and t2.

% Some dummy values; delete these once you start writing your own code.
rg1 = [Re + 450*1000; Re + 550*1000; Re + 650*1000];
vg1 = [-2000; 1000; 4000];

end