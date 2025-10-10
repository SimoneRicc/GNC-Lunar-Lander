function I_total = total_inertia(masses, centers_of_mass, inertias)
% masses:            1xN (o Nx1) masse dei corpi
% centers_of_mass:   3xN, ogni colonna è il CoM del corpo nel frame globale
% inertias:          1xN cell, ogni cella una 3x3 = tensore d’inerzia centroidale (stesso frame)

    masses = masses(:);              % Nx1
    total_mass = sum(masses);        

    % Centro di massa complessivo (colonne*masse)
    center_of_mass_total = (centers_of_mass * masses) / total_mass;  % 3x1

    I_total = zeros(3,3);
    N = numel(masses);

    for i = 1:N
        Ic = inertias{i};            % 3x3 centroidale nel frame globale
        r  = centers_of_mass(:,i) - center_of_mass_total;     % 3x1
        I_total = I_total + Ic + masses(i) * ( (r.'*r)*eye(3) - (r*r.') );
    end

    % Simmetrizzazione numerica
    I_total = 0.5*(I_total + I_total.');
end
