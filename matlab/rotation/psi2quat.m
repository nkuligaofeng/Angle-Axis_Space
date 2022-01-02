function quat = psi2quat(psi)
    angle_in_radian = norm(psi);
    if angle_in_radian == 0
        w = 1; x = 0; y = 0; z = 0;
    else
        w = cos(angle_in_radian / 2.0);
        x = sin(angle_in_radian / 2.0) * psi(1) / angle_in_radian;
        y = sin(angle_in_radian / 2.0) * psi(2) / angle_in_radian;
        z = sin(angle_in_radian / 2.0) * psi(3) / angle_in_radian;
    end
    quat = [w; x; y; z];
end