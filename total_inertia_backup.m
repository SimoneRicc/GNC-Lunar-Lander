function I_total = total_inertia(masses, centers_of_mass, inertias)
    % masses: vettore delle masse dei corpi
    % centers_of_mass: matrice 3xN dei centri di massa
    % inertias: cell array delle matrici di inerzia dei corpi

    % Calcola il centro di massa totale
    total_mass = sum(masses);
    center_of_mass_total = sum(masses .* centers_of_mass, 2) / total_mass;

    % Inizializza la matrice di inerzia totale
    I_total = zeros(3, 3);

    for i = 1:length(masses)
        % Trasforma la matrice di inerzia al centro di massa totale
        r = centers_of_mass(:, i) - center_of_mass_total; % vettore di traslazione
        I_transformed = inertias(i) + masses(i) * (norm(r)^2 * eye(3) - r * r'); % Teorema degli assi paralleli
        I_total = I_total + I_transformed;
    end
end