clear;
% case1: Numerical calculate 
% define the variables
syms Alr pp pm su s1 s2 s3 s4 dp p N1 N2 N3 N4 DP
% Define Alr as a function under specific conditions
f_alr = -(sqrt((s1 - s2 - s3 + s4) * (s1 + s2 + s3 + s4) * ...
    (s1^2 - s2^2 + 2 * (1 - 2 * dp^2) * s2 * s3 - s3^2 + ...
    2 * (-1 + 2 * dp^2) * s1 * s4 + s4^2)) / ...
    ((s1 - s2 - s3 + s4) * (s1 + s2 + s3 + s4)));
% Find the partial derivative of Alr
grad_alr = gradient(f_alr, [s1, s2, s3, s4, dp]);
% Calculate the delta Alr value 
Pro_delta_Alr = sqrt( ...
     (grad_alr(1) .* (1 ./ sqrt(N1)) .* su .* (1 - (pp+dp) * (pm+dp) + Alr * (pp - pm))).^2 ...
     + (grad_alr(2) .* (1 ./ sqrt(N2)) .* su .* (1 + (pp - dp) * (pm+dp) + Alr * ((-pp + dp) - (pm+dp)))).^2 ...
     + (grad_alr(3) .* (1 ./ sqrt(N3)) .* su .* (1 + (pp+dp) * (pm - dp) + Alr * ((pp+dp) + pm - dp))).^2 ...
     + (grad_alr(4) .* (1 ./ sqrt(N4)) .* su .* (1 - (pp - dp) * (pm - dp) + Alr * (-pp + pm))).^2 ...
     + (grad_alr(5) .* (DP).*pm).^2 ...
    );
% Symbolic substitution calculation
Pro_delta_Alr1 = subs(Pro_delta_Alr, {'s1', 's2', 's3', 's4'}, ...
    {su .* (1 - (pp+dp) * (pm + dp) + Alr * (pp - pm)), ...
    su .* (1 + (pp - dp) * (pm + dp) + Alr * ((-pp + dp) - (pm+dp))), ...
    su .* (1 + (pp + dp) * (pm - dp) + Alr * ((pp+dp) + pm - dp)), ...
    su .* (1 - (pp - dp) * (pm - dp) + Alr * (-pp + pm))}); 
% Define numeric values
Alr_val = 0.16;
pp_val = 0.4;
pm_val = 0.7;
su_val = 2;
N1_val = 0.25 * 10^9;
N2_val = 0.25 * 10^9;
N3_val = 0.25 * 10^9;
N4_val = 0.25 * 10^9;

% Define ranges for dp_val and DP_val
dp_values = linspace(0.01, 0.1, 10); 
DP_values = logspace(-3, -1, 10);   

% Initialize result matrix
Pro_delta_Alr_results = zeros(length(dp_values), length(DP_values));

% Perform numerical calculations in loops
for i = 1:length(dp_values)
    for j = 1:length(DP_values)
        dp_val = dp_values(i);
        DP_val = DP_values(j);
        
        % Substitute the numeric values into Pro_delta_Alr1
        Pro_delta_Alr_num = subs(Pro_delta_Alr1, ...
            {'su', 'pp', 'pm', 'Alr', 'dp', 'N1', 'N2', 'N3', 'N4', 'DP'}, ...
            {su_val, pp_val, pm_val, Alr_val, dp_val, N1_val, N2_val, N3_val, N4_val, DP_val});
        
        % Store the numerical result
        Pro_delta_Alr_results(i, j) = double(Pro_delta_Alr_num);
    end
end

% Generate 3D plot with exponential colorbar
[DP_mesh, dp_mesh] = meshgrid(DP_values, dp_values);
figure;
surf(dp_mesh, log10(DP_mesh), Pro_delta_Alr_results);
xlabel('$\Delta P$', 'Interpreter', 'latex');
ylabel('$lg\frac{\delta P}{P}$', 'Interpreter', 'latex');
zlabel('$\delta A_{LR}$', 'Interpreter', 'latex');
title('The $\delta A_{LR}$ get from propagation of error', 'Interpreter', 'latex')

% Use logarithmic color scaling
c = colorbar;
set(c, 'Ticks', logspace(log10(min(Pro_delta_Alr_results(:))), log10(max(Pro_delta_Alr_results(:))), 5), ...
       'TickLabels', arrayfun(@(x) sprintf('10^{%.1f}', log10(x)), ...
       logspace(log10(min(Pro_delta_Alr_results(:))), log10(max(Pro_delta_Alr_results(:))), 5), 'UniformOutput', false));

% Add shading and color map for better gradient visualization
shading interp;
colormap jet;

% Add gradient plot (2D contour map)
figure;
contourf(dp_mesh, log10(DP_mesh), Pro_delta_Alr_results, 20, 'LineColor', 'none');
xlabel('$\Delta P$', 'Interpreter', 'latex');
ylabel('$lg\frac{\delta P}{P}$', 'Interpreter', 'latex');
title('The $\delta A_{LR}$ get from propagation of error', 'Interpreter', 'latex');
colorbar;
colormap jet;


