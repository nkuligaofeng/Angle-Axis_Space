function quat = rotm2quat(R)
    psi = LogRotation(R);
    q = psi2quat(psi);
    quat.s = q(1);
    quat.v = q(2:4);
end