%% The global warming-cooling effect numerical modelling from the obtained oscillations 
%% using the non-linear logistic delay-differential equation (LODE) algorithm  
% Solve equation system (1,2) from Rial 2004 (modified)
%% Any use of this code must refer to
% Kravchinsky,V.A., Zhang,R., Borowiecki,R., Goguitchaichvili,A., Czarnecki,J., Czarnecki,A., Boers,N., Berger,A., & van der Baan,M. (2025)
% Millennial Cycles in Greenland and Antarctic Ice Core Records: Evidence of Astronomical Influence on Global Climate
% Journal of Geophysical Research: Atmospheres
% http//
% Questions/suggestions to vadim@ualberta.ca

clc; clear; close all;

%% Constants
time_offset = -9;   % If time_offset more negative, peaks shifted to more remote past
heat_capacity = 16.7 / 1000;
initial_temperature = 10; % Degrees Celsius
solar_cont_quarter = 1368 / 4;
initial_ice = 0.9;
tmin = -160;
tmax = 50; % Unit: kiloyears
dt = 0.005;
Npoints = floor((tmax - tmin) / dt);
% Albedo Maximum Value
Albedo_Max = 0.65;

if time_offset <= 0
    tmin = tmin + time_offset;
    nmin = floor(-time_offset / dt);
    nmax = Npoints;
else
    tmax = tmax + time_offset;
    nmin = 1;
    nmax = Npoints + floor(-time_offset / dt);
end

time_points = linspace(tmin, tmax, Npoints);
ice = zeros(size(time_points));
ice(1) = initial_ice;
temperature = zeros(size(time_points));
temperature(1) = initial_temperature;


%% Define functions
function val = clouds(L)
    Lminus = 0.9; Lplus = 1.5;
    if L >= Lplus
        val = 1;
    elseif L < Lminus
        val = 0;
    else
        val = 0.5; % Lminus <= L < Lplus
    end
end

function val = Albedo(L, Albedo_Max)
    val = Albedo_Max * tanh(L / 2);
end

function val = temperature_RHS(temperature, ice, n, solar_cont_quarter, Budyko_a, Budyko_b, Budyko_a1, Budyko_b1, heat_capacity, Albedo_Max)
    val = (solar_cont_quarter * (1 - Albedo(ice(n-1), Albedo_Max)) - Budyko_a - Budyko_b * temperature(n-1) + ...
        (Budyko_a1 + Budyko_b1 * temperature(n-1)) * clouds(ice(n-1))) / heat_capacity;
end

function val = forcing(t)
    periods = [22, 11.4, 5.5, 2.75, 1.4, 1];
    amplits = [0.0, 0.22, 0.2, 0.12, 0.1, 0.1];
    val = sum(amplits .* cos(2 * pi * t ./ periods));
end

function val = K(n, dt, temperature, AMPLITUDE)
    time = n * dt;
    val = 1 + AMPLITUDE * (1 + forcing(time)) * temperature(n-1);
end

function val = RHS(temperature, L, n, delay_n, time_constant_mu, both, dt, AMPLITUDE)
    point_of_influence = n - 1 - delay_n;
    if point_of_influence >= 1
        if both
            val = time_constant_mu * L(point_of_influence) * (1 - L(point_of_influence) / K(n, dt, temperature, AMPLITUDE));
        else
            val = time_constant_mu * L(n-1) * (1 - L(point_of_influence) / K(n, dt, temperature, AMPLITUDE));
        end
    else
        val = time_constant_mu * L(n-1) * (1 - L(n-1) / K(n, dt, temperature, AMPLITUDE));
    end
end

%% Parameters
scaling = 30; % Factor by which Greenland ice is smaller than total, after Rial:2004
time_constant_mu = scaling * 1/14; % 1/.48 %1./11.
Budyko_a = 223;
Budyko_b = 2.2;
Budyko_a1 = 47.8;
Budyko_b1 = 1.6;
AMPLITUDE = 0.05;
BOTH_DELAYS_Q = true;
ABOVE_BIFFUR = 1.1;

delay_tau = ABOVE_BIFFUR * (pi / 2) / time_constant_mu;
delay_n = floor(delay_tau / dt);

%% Run the model
for n = 2:length(ice)
    temperature(n) = temperature(n-1) + dt * temperature_RHS(temperature, ice, n, solar_cont_quarter, Budyko_a, Budyko_b, Budyko_a1, Budyko_b1, heat_capacity, Albedo_Max);
    ice(n) = ice(n-1) + dt * RHS(temperature, ice, n, delay_n, time_constant_mu, BOTH_DELAYS_Q, dt, AMPLITUDE);
end

%% Ensure `nmin` and `nmax` are within valid index range
nmin = max(1, floor(nmin));
nmax = min(length(time_points), floor(nmax));

% temperature = temperature - min(temperature);
% ice = ice-min(ice);

%% Plot results
figure (1)
set(gcf, 'Position',  [0, 0, 1000, 600])

subplot(2,1,1);
plot(-time_points(nmin:nmax), ice(nmin:nmax), 'b', 'LineWidth', 1.5);
ylabel('Ice');
% ylim([1, 2]);
axis([-20 150 1.3 1.45]); 
xticks(-20:10:150); 
grid on;

subplot(2,1,2);
plot(-time_points(nmin:nmax), temperature(nmin:nmax), 'r', 'LineWidth', 1.5);
xlabel('ky before present');
ylabel('T [Celsius]');
% ylim([0, 40]);
axis([-20 150 4 12]); 
xticks(-20:10:150); 
grid on;

%% Requires R2020a or later
exportgraphics(gcf,'vectorfig.pdf','ContentType','vector')

%% Save to CSV
csv_file = "timeBP_temperature_ice.csv";
data = [-time_points(nmin:nmax)', temperature(nmin:nmax)', ice(nmin:nmax)'];
writematrix(data, csv_file);
disp("The data export to .csv is done.");
