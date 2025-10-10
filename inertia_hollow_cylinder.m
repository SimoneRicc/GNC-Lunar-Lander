function I = inertia_hollow_cylinder(r, R, h)
    Ixx = (1/12) * (3*(r^2+R^2) + h^2);
    Iyy = Ixx;
    Izz = 0.5 * (r^2+R^2);
    I = diag([Ixx, Iyy, Izz]);
end