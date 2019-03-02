function C1_matrix = C1(th)
% This is a ``1'' axis DCM.

C1_matrix = [1     0        0;
             0  cos(th)  sin(th);
             0 -sin(th) cos(th)  ];