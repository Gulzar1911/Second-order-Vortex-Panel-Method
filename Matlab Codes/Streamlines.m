



%q6_assigned

    % Parameters
    c = 1.0;                   % Chord length
    epsilon_over_c = 0.11;     % Maximum camber to chord ratio
    angle_of_attack = 7;       % Angle of attack in degrees
    m = 100;                   % Number of panels
    Vinf=1;

    filename = sprintf('airfoil_data_100_0.11.txt');
    generate_and_save_airfoil_data(m, epsilon_over_c, c, filename);
    
    
    % Load the airfoil data
    data = load("airfoil_data_100_0.11.txt", '-ascii'); % Use '-ascii' to read text files
    Xb = data(:, 1);
    Yb = data(:, 2);

    [c_l, Cp , s] = VORTEX_PANEL_FUNCTION(X, Y, Vf, aoa);
  


 

function generate_and_save_airfoil_data(m, epsilon_over_c, c, filename)
    % Ensure m is even
    if mod(m, 2) ~= 0
        error('m must be an even number');
    end
    
    % Number of points on each surface
    num_points_top = m / 2 + 1;
    num_points_bottom = m / 2;
    
    % Initialize arrays for x and y coordinates
    X_top = zeros(1, num_points_top);
    Y_top = zeros(1, num_points_top);
    X_bottom = zeros(1, num_points_bottom);
    Y_bottom = zeros(1, num_points_bottom);
    
    % Calculate top surface coordinates
    for j = 1:num_points_top
        % Calculate normalized x-coordinate for top surface
        X_j_over_c = 2 * (j - 1) / m;
        
        % Calculate actual x-coordinate for top surface
        X_top(j) = X_j_over_c * c;
        
        % Calculate y-coordinate using the given formula for top surface
        Y_top(j) = 4 * epsilon_over_c * X_j_over_c * (1 - X_j_over_c);
    end
    
    % Calculate bottom surface coordinates
    for j = (m / 2 + 2):(m + 1)
        % Calculate normalized x-coordinate for bottom surface
        X_j_over_c = 2 * (m - j + 1) / m;
        
        % Calculate actual x-coordinate for bottom surface
        X_bottom(j - (m / 2 + 1)) = X_j_over_c * c;
        
        % y-coordinate is zero for the bottom surface
        Y_bottom(j - (m / 2 + 1)) = 0;
    end
    
    % Combine top and bottom surface coordinates
    Xb = [X_top, X_bottom];
    Yb = [Y_top, Y_bottom];
   
    % Save data to file
    data = [Xb', Yb'];
    save(filename, 'data', '-ascii');
end



function [c_l, s, Cp] = VORTEX_PANEL_FUNCTION(X, Y, Vf, aoa)
    %% Variables declarations
    Mp = length(X); % Number of boundary points
    M = Mp - 1; % Number of panels

    %% For Coordinates (x,y) of control point and panel length S computed for each of the vortex panels
    x = zeros(1, M);
    y = zeros(1, M);
    s = zeros(1, M);
    phi = zeros(1, M);
    eta = zeros(1, M);
    beta = zeros(1, M);
    sine = zeros(1, M);
    cosine = zeros(1, M);
    RHS = zeros(1, M);

    for j = 1:M
        x(j) = 0.5 * (X(j) + X(j+1));
        y(j) = 0.5 * (Y(j) + Y(j+1));
        s(j) = sqrt((X(j+1) - X(j))^2 + (Y(j+1) - Y(j))^2);

        phi(j) = atan2((Y(j+1) - Y(j)), (X(j+1) - X(j))); % Calculate phi
        eta(j) = phi(j) + (pi/2);
        beta(j) = eta(j) - deg2rad(aoa); % Convert aoa to radians
        sine(j) = sin(phi(j)); % Calculate sine(phi)
        cosine(j) = cos(phi(j)); % Calculate cosine(phi)
        RHS(j) = 2 * pi * Vf * cos(beta(j)); % RHS represents the right-hand side
    end

    K1 = zeros(M, M);
    K2 = zeros(M, M);

    for i = 1:M
        for j = 1:M
            % Calculate the geometric constants
            A = -(x(i) - X(j)) * cosine(j) - (y(i) - Y(j)) * sine(j);
            B = (x(i) - X(j))^2 + (y(i) - Y(j))^2;
            C = -cos(phi(i) - phi(j));
            D = (x(i) - X(j)) * cosine(i) - (y(i) - Y(j)) * sine(i);
            E = sqrt(abs(B - A^2));
            F = atan(A / E);
            G = log((s(j)^2 + 2 * A * s(j) + B) / B);
            H = atan2(s(j) + A, E);
            I = D - 2 * A * C;

            if E > 0
                K1(i, j) = C / 2 * G + (D - A * C) / E * (H - F);
                K2(i, j) = C + (I / (2 * s(j))) * G - ((B * C + A * I) / (E * s(j))) * (H - F);
            elseif E == 0
                K1(i, j) = C / 2 * G + (D - A * C) * ((1 / (s(j) + A)) - (1 / A));
                K2(i, j) = C + (I / (2 * s(j))) * G - ((B * C + A * I) / s(j)) * ((1 / (s(j) + A)) - (1 / A));
            end
        end
    end

    %% Compute Influence Coefficients
    AN = zeros(M, M);

    for i = 1:M % Excluding Mth panel Boundary Condition
        AN(i, 1) = K2(i, M) - K2(i, 1) + K1(i, 1);
        for j = 2:M-1
            AN(i, j) = K1(i, j-1) - K2(i, j) + K1(i, j);
        end
    end

    AN(M, :) = 1;

    %% Solve the Linear System
    Gamma = AN \ RHS';

    %% Calculate the Sectional Coefficient of Lift
    c = zeros(1, M);
    for j = 1:M
        if j == 1
            c(j) = (Gamma(j+1) + Gamma(j)) / 2 * s(j);
        else
            c(j) = (Gamma(j) + Gamma(j-1)) / 2 * s(j);
        end
    end
    
    L = sum(c) * Vf;
    c_l = L / (0.5 * Vf^2 * sum(c));

    %% Calculate pressure Coefficient
    Cp = 1 - (Gamma / Vf).^2;
    
    %% Plot the streamlines
    [x_1, y_1] = meshgrid(-0.5:0.002:1.5, -0.5:0.001:0.5);
    psi = zeros(1001, 1001);
    for i = 1:length(x_1)
        for j = 1:length(y_1)
             psi(i, j) = Vf .* y_1(i, j);
            for k = 1:(length(Gamma) - 1)
                psi(i, j) =  psi(i, j) + (Gamma(k) ./ (2 * pi)) .* (log((y_1(i, j) - y(k)).^2 + (x_1(i, j) - x(k)).^2));
            end
        end
    end
    figure;
    contour(x_1, y_1, psi, 250); % contour code to plot 250 streamlines
    hold on
    plot(X, Y) % Plot airfoil based on the data
    set(gca,'FontWeight','bold')
    set(gca,'FontSize', 10)
    hold off
    aoa = (aoa * 180) / pi; % Convert back to degrees
    name_1 = strcat('Streamlines for \alpha= ', num2str(aoa));
    title(name_1)
    xlabel ('x-axis')
    ylabel ('y-axis')
    grid on
    axis image
    
    
end
