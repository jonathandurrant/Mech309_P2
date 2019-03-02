function C2_matrix = C2(th)
% This is a ``2'' axis DCM.

C2_matrix = [  cos(th)  0 -sin(th);
                  0     1     0   ;
               sin(th)  0 cos(th) ];