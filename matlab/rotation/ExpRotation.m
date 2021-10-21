%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %   get a 3x3 rotation matrix from a rotation vector %

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input:   w ----- 3x1 rotation vector
%  Output:  R ------ rotation matrix
%  Made by: Gaofeng Li (PhD student at Nankai University)
%  Date:    Nov 17, 2016
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = ExpRotation(w)
    [m,n]=size(w);
    if m ~= 3 || n ~= 1
        disp('w is not a 3x1 vector');                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    end
    w_norm = sqrt( w(1)^2 + w(2)^2 + w(3)^2 );
    if w_norm ~= 0
        skewMatrix = [0, -w(3)/w_norm, w(2)/w_norm;
                  w(3)/w_norm, 0, -w(1)/w_norm;
                  -w(2)/w_norm, w(1)/w_norm, 0];
        R = eye(3) + skewMatrix * sin(w_norm) + (1-cos(w_norm))*(skewMatrix*skewMatrix);
    else
        R = eye(3);
    end