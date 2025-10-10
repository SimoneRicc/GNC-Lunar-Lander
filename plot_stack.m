function plot_stack(r_LLM,h_LLM,r_ST,h_ST,H_geo,W_geo,L_geo, LLM_cog, ST_cog, GEOSAT_cog, r_cog, comp, name_title)
    % Parametri geometrici
    radius_LLM = r_LLM;
    height_LLM = h_LLM;

    radius_ST = r_ST;
    height_ST = h_ST;

    H = H_geo;  % satellite height (Z)
    W = W_geo;
    L = L_geo;

    z_LLM = 0;
    z_ST  = z_LLM + height_LLM;
    z_SAT = z_ST  + height_ST;

    figure;
    hold on; 
    grid on; 
    axis equal;
    view(3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(name_title);

    % Draw elements
    h1 = draw_cylinder(radius_LLM, height_LLM, [0 0 z_LLM], [1 0 0]);
    h2 = draw_cylinder(radius_ST, height_ST, [0 0 z_ST-comp], [0 1 0]);
    h3 = draw_box(L, W, H, [ -L/2, -W/2, z_SAT-comp], [0 0 1]);
    h4 = draw_point(LLM_cog,0.2, [1 0 0]);
    h5 = draw_point(ST_cog,0.2, [0 1 0]);
    h6 = draw_point(GEOSAT_cog,0.2, [0 0 1]);
    h7 = draw_point(r_cog, 0.6, [0 0 0]);

    legend([h1, h2, h3, h4, h5 ,h6 ,h7], {'LLM ','ST','GEOSAT', 'LLM CoG', 'ST CoG','GEOSAT CoG','Multi-Body CoM'});
end

% === Auxiliary functions ===

function h = draw_cylinder(r, h_cyl, origin, color)
    % Superficie laterale
    N = 30;
    [X, Y, Z] = cylinder(r, N);
    Z = Z * h_cyl;
    
    X = X + origin(1);
    Y = Y + origin(2);
    Z = Z + origin(3);
    
    % Superficie laterale
    h = surf(X, Y, Z, ...
        'FaceAlpha', 0.1, 'FaceColor', color, 'EdgeColor', 'none');

    % Tappi superiore e inferiore
    theta = linspace(0, 2*pi, N);
    x_circle = r * cos(theta) + origin(1);
    y_circle = r * sin(theta) + origin(2);

    fill3(x_circle, y_circle, origin(3)*ones(1,N), color, ...
        'FaceAlpha', 0.1, 'EdgeColor', 'k');
    fill3(x_circle, y_circle, (origin(3)+h_cyl)*ones(1,N), color, ...
        'FaceAlpha', 0.1, 'EdgeColor', 'k');
end

function h = draw_box(L, W, H, origin, color)
    [X, Y, Z] = ndgrid([0 L], [0 W], [0 H]);
    verts = [X(:), Y(:), Z(:)];
    verts = verts + origin;
    K = convhull(verts(:,1), verts(:,2), verts(:,3));
    h = trisurf(K, verts(:,1), verts(:,2), verts(:,3), ...
        'FaceColor', color, 'FaceAlpha', 0.1, 'EdgeColor', 'k');
    % === PANNELLI SOLARI CENTRATI, ROTATI E TRASLATI ===
    panel_W = 3;     % lato corto (altezza su Z)
    panel_L = 6;     % lato lungo (esteso su Y)
    panel_t = 0.05;    % spessore visivo (trascurabile)
    panel_offset_Y = 0.1;  % distanza extra dal box per non intersecare

    panel_color = color; 

    % Centro del box
    x_center = origin(1) + L/2;
    y_left   = origin(2) - panel_offset_Y;          % sposta più indietro
    y_right  = origin(2) + W + panel_offset_Y;      % sposta più avanti
    z_center = origin(3) + H/2;

    % Coordinate geometriche pannelli
    y1 = -panel_L/2-10;
    y2 = +panel_L/2+10;
    z1 = z_center - panel_W/2;
    z2 = z_center + panel_W/2;

    % === Pannello sinistro (faccia -Y, spostato all’indietro) ===
    verts_L = [
        x_center, y_left + y1, z1;
        x_center, y_left + y2, z1;
        x_center, y_left + y2, z2;
        x_center, y_left + y1, z2
    ];
    patch('Vertices', verts_L, 'Faces', [1 2 3 4], ...
          'FaceColor', panel_color, 'FaceAlpha', 0.1, 'EdgeColor', 'k');

    % === Pannello destro (faccia +Y, spostato in avanti) ===
    verts_R = verts_L;
    verts_R(:,2) = verts_R(:,2) + (y_right - y_left);  % traslazione Y
    patch('Vertices', verts_R, 'Faces', [1 2 3 4], ...
          'FaceColor', panel_color, 'FaceAlpha', 0.1, 'EdgeColor', 'k');

end

function h = draw_point(center, radius, color)
    % Disegna una piccola sfera 3D
    [X, Y, Z] = sphere(15);  % risoluzione della mesh
    X = X * radius + center(1);
    Y = Y * radius + center(2);
    Z = Z * radius + center(3);

    h = surf(X, Y, Z, ...
        'FaceColor', color, ...
        'FaceAlpha', 1.0, ...
        'EdgeColor', 'none');
end