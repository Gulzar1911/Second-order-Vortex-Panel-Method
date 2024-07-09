

% Parameters
c = 1.0;                   % Chord length
epsilon_over_c = 0.11;     % Maximum camber to chord ratio
aoa = 7;  % Angle of attack in degrees
m = 100;  % Number of panels
Vinf = 1;

  
% Generate airfoil data

    filename = sprintf('airfoil_data_100_0.11.txt');
    generate_and_save_airfoil_data(m, epsilon_over_c, c, filename);
    
    
    % Load the airfoil data
    data = load("airfoil_data_100_0.11.txt", '-ascii'); % Use '-ascii' to read text files
    Xb = data(:, 1);
    Yb = data(:, 2);
    
  
    [c_l, Cp, s] = VORTEX_PANEL_FUNCTION(X, Y, c, aoa);





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
    c_l = L / (0.5 * Vf^2 * sum(s));

    %% Calculate pressure Coefficient
    Cp = 1 - (Gamma / Vf).^2;

    % Separate CP for upper and lower surfaces
    Cp_upper = Cp(1:M/2);    % Upper surface Cp values
    Cp_lower = Cp(M/2+1:end);  % Lower surface Cp values

    X_upper = x(1:M/2);
    X_lower = x(M/2+1:end);

    %% Plot coefficient of pressure graph
    figure;
    plot(X_upper ./ c(1), Cp_upper, 'b-', 'LineWidth', 1, 'DisplayName', 'Upper Surface'); % Plot pressure coefficient distribution
    hold on;
    plot(X_lower ./ c(1), Cp_lower, 'r-', 'LineWidth', 1, 'DisplayName', 'Lower Surface'); % Plot pressure coefficient distribution
    title(['Pressure Coefficient vs x/c for \alpha = ', num2str(aoa), ' degrees']);
    xlabel('x/c');
    ylabel('C_p');
    set(gca, 'YDir', 'reverse'); % Flip the y-axis for Cp values
    set(gca, 'FontWeight', 'bold', 'FontSize', 10);
    legend show;
    grid on;
end
