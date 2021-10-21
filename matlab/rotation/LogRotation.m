%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %   get a 3*1 angle-axis vector from a Rotation Matrix %

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input:   R ----- Rotation Matrix
%  Output:  w ------ The 3*1 angle-axis vector. Its L2-norm is constrainted less than Pi;
%  Made by: Gaofeng Li (PhD student at Nankai University)
%  Date:    Feb 22, 2017
%           Modified in Oct. 10, 2017. the special case 1 for "R is I3" is
%           changed from R == eye(3) to norm(R - eye(3),'fro') < 0.0001
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = LogRotation(R)
    [m,n]=size(R);
    if m ~= 3 || n ~= 3
        disp('R is not a 3x3 matrix, please check again');
    else
        if norm(R'*R - eye(3),'fro') > 0.00001
            disp('R^T*R is not a identity matrix, please check again');
        else
            if abs(det(R)-1) > 0.00001
                disp('det(R) is not 1, please chech again');
            else
              %% %%%%%%%%%%%%%%%%%The core of this algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
                if norm(R - eye(3),'fro') < 0.0001 %special case 1:the rotation matrix is I3
                    w = [0; 0; 0;];
                else
                    Omega = 0.5*(R - R'); %skew-symmetric matrix * 0.5
                    wTemp = [-Omega(2,3); Omega(1,3); -Omega(1,2)];
                    if norm(wTemp) < 0.00001 %special case 2: the angle is pi;(sometimes it may not be exactly 0)
                         %根据R=I+sin(pi)*wUnit^ + (1-cos(pi))/(pi^2)*wUnit^*wUnit^;
                         %wUnit是w的单位向量对应的反对称矩阵；因为w和-w是表示同一个矩阵。所以我们总是取wx为正的那一个
                         wx = sqrt( (R(1,1)+1)*pi^2/2 );
                         if wx == 0
                            wy = sqrt( (R(2,2)+1)*pi^2/2 );
                            if wy == 0
                                wz = sqrt( (R(3,3)+1)*pi^2/2 );
                            else
                                wz = R(2,3)*pi^2/(2*wy);
                            end
                         else %wx不为0的情况
                             wy = R(1,2)*pi^2/(2*wx);
                             wz = R(1,3)*pi^2/(2*wx);
                         end
                         w = [wx; wy; wz];
                    else
                        wNorm = acos( ( R(1,1)+R(2,2)+R(3,3)-1 )/2 );
                        w = (wNorm/norm(wTemp))*wTemp;
                    end
                end
              %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    end