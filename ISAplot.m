% ISA plot script - Gets properties for a wanted set of altitudes and plots them  
%
% See also
%   ISAfunction.m
%
% References
%   <Reference List like Papers>
%   ...
%
% Authors
%   Niklas Gierse           [nikgi434]
%   Leonhard Muehlstrasser  [leomu719]
%
% License
%   This program is part of an academic exercise for the course TMAL02,  
%   Linköping University, year 2023. The program/template is therefore   
%   free for non-commercial academic use.  
%
% Code History
%   https://github.com/ngiersetum/tmal02_lab1

%% The executable code starts here:
% Calculate atmospheric properties with ISAfunction.m and plot them

% Altitudes we want values for. Steps of 50m produce nice plots
altitudes = 0:50:86000;
    
% Calculate
[T, p, rho, a, mu, id1, id2] = ISAfunction(altitudes);

% Scale all the values against the sea level value to plot them all
% together

% Same values as in ISAfunction.m
T0 = 288.15; %[K]
p0 = 101325; % [Pa]
rho0 = 1.2250; %[kg/m^3]
a0 = 340.294; %[m/s]
mu0 = 1.7894e-5; %[kg/m*s]

TScaled = T ./ T0;
pScaled = p ./ p0;
rhoScaled = rho ./ rho0;
aScaled = a ./ a0;
muScaled = mu ./ mu0;

% Plot it all
figure
hold on
grid on
xlim([0 1])

plot(TScaled(:), altitudes(:), 'Color', 'red', 'LineWidth', 1)
plot(pScaled(:), altitudes(:), 'b')
plot(rhoScaled(:), altitudes(:))
plot(aScaled(:), altitudes(:))
plot(muScaled(:), altitudes(:))

% Styling
xlabel('T/T_{0} , P/P_{0} , \rho/\rho_{0} , a/a_{0} , \mu/\mu_{0}');
ylabel('Z [m]');
title('ISA Atmospheric properties');
legend('T/T_{0}', 'P/P_{0}', '\rho/\rho_{0}', 'a/a_{0}', '\mu/\mu_{0}', 'Location', 'NorthWest');

