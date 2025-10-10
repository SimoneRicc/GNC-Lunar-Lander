function I = inertia_box(a, b, c)
    % Momento d'inerzia di un parallelepipedo (a, b, c: dimensioni)
    Ixx = (1/12) * (b^2 + c^2);
    Iyy = (1/12) * (a^2 + c^2);
    Izz = (1/12) * (a^2 + b^2);
    I = diag([Ixx, Iyy, Izz]);
end