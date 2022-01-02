function R = quat2rotm(quat)
    %quat: w, x, y, z
    [m, n] = size(quat);
    if ~( (m == 4 && n == 1) || (m == 1 && n == 4) )
        disp('The size of quaternion is incorrect: [m,n] = [', num2str(m), ', ', num2str(n), ']');
    end
    q_norm = sqrt( quat(1)^2 + quat(2)^2 + quat(3)^2 +quat(4)^2);
    if abs( q_norm - 1 ) > 0.0001
        disp('Warning!!!!---The input quaternion is not a unit quaternion.---');
        quat = quat / q_norm;
        disp('It has been normalized');
    end
    r11 = 1 - 2 * quat(3)^2 - 2 * quat(4)^2;
    r12 = 2 * quat(2) * quat(3) - 2 * quat(1) * quat(4);
    r13 = 2 * quat(2) * quat(4) + 2 * quat(1) * quat(3);
    r21 = 2 * quat(2) * quat(3) + 2 * quat(1) * quat(4);
    r22 = 1 - 2 * quat(2)^2 - 2 * quat(4)^2;
    r23 = 2 * quat(3) * quat(4) - 2 * quat(1) * quat(2);
    r31 = 2 * quat(2) * quat(4) - 2 * quat(1) * quat(3);
    r32 = 2 * quat(3) * quat(4) + 2 * quat(1) * quat(2);
    r33 = 1 - 2 * quat(2)^2 - 2 * quat(3)^2;
    R = [r11 r12 r13;
        r21 r22 r23;
        r31 r32 r33];
    
end