%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %   return the geodesic distance bewteen two rotation matrixs  %

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input:   R1, R2 ----- the two Rotation Matrix
%  Output:  d ------ The geodesic distance of two rotation matrix
%  Made by: Gaofeng Li (PhD student at Nankai University)
%  Date:    Feb 24, 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w, d] = DistanceRotations(R1,R2)
    w = LogRotation(R1'*R2);
    d = norm(w);