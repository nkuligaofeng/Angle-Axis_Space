%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %   get a the foot of perpendicular %

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input:   R_tar ----- target orientation
%           R_cur ----- current orientation
%           omega ----- the tool direction
%  Output:  R_perp ------ the foot of perpendicular
%  Made by: Gaofeng Li (PhD student at Nankai University)
%  Date:    April 25, 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R_perp = GetRPerp(R_tar, R_cur, omega)
omega_1 = LogRotation(R_tar'*R_cur);
omega_2 = 2*(omega'*omega_1)*omega - omega_1;
omega_bar = LogRotation(R_cur'*R_tar*ExpRotation(omega_2));
if norm(omega_bar)<0.0001
    R_perp = R_cur;
else
    R_c = R_tar'*R_cur*ExpRotation(0.5*omega_bar);
    if norm(cross(LogRotation(R_c),omega)) < 0.0000001
        R_perp = R_tar*R_c;
    else
        if norm(omega_bar) ~=0
            R_perp = R_cur*ExpRotation( -( ( 2*pi-norm(omega_bar) )/( 2*norm(omega_bar) ) )*omega_bar );
        end
    end
end