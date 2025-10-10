function I = inertia_cylinder(r, h)
    Ixx = (1/12) * (3*r^2 + h^2);
    Iyy = Ixx;
    Izz = 0.5 * r^2;
    I = diag([Ixx, Iyy, Izz]);
end