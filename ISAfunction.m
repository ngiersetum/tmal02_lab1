function [T, P, rho, a, mu, liuID1, liuID2] = ISAfunction(Z)
% ISA function - Calculate atmospheric properties for one or several
% altitudes
%
% Synopsis
%   [T P rho a mu liuid1 liuid2] = ISAfunction(Z)
%
% Inputs
%   (scalar or vector) Z      Geometric altitude in [m]
%
% Outputs
%   (scalar or vector) T      Temperature, K
%   (scalar or vector) P      Pressure, Pa
%   (scalar or vector) rho    Density, kg/m^3
%   (scalar or vector) a      Sound speed, m/s
%   (scalar or vector) mu     Dynamic viscosity, kg/(m s)
%
% References
%   <Reference List like Papers>
%   ...
%
% Authors
    liuID1 = "nikgi434"; % author 1 name / <LiU ID 1>
    liuID2 = "leomu719"; % author 2 name / <LiU ID 2>
%
% License
%   This program is part of an academic exercise for the course TMAL02,
%   Linköping University, year 2023. The program is therefore free for 
%   non-commercial academic use.
%
% Code History
%   https://github.com/ngiersetum/tmal02_lab1

%% The executable code starts here:
% Define constants and reference values
    
    g0 = 9.80665; %[m/s²]
    r0 = 6356766; %[m]
    R = 8.31432e3; %[(N*m)/(kmol*K)]
    gamma = 1.4; %dimensionless
    beta = 1.458e-6; %[kg/(s*m*K^1/2)]
    bs = [0, 1, 2, 3, 4, 5, 6, 7];
    S = 110; %[K] (Sutherland Constant)
    Hb =[0, 11000, 20000, 32000, 47000, 51000, 71000, 84852]; %[m']
    bet = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.0020]; %[K/m]
    Tb =[288.15, 216.650, 216.650, 228.650, 270.650, 270.560, 214.650, 186.8673]; %[K]
    pb = [101325, 22632, 5474.8, 868.01, 110.9, 66.938, 3.9564, 0.39814]; % [Pa]
    rhob = [1.2250]; %[kg/m³]
    M0 = 28.9644; %[kg/kmol]
    ab = [340.294]; %[m/s]
    mub = [1.7894e-5]; %[kg/ms]
    
    
% Determine number of input altitudes
    nInput = length(Z);
    
% Convert input (geometric altitude) to geopotential altitude for
% calculations. Conversion is done as an individual operation for every
% altitude in the loop later
    g = @(Z) g0 .* (r0./(r0+Z)).^2;   % function for calculation of g

% Initialize arrays just so they don't change size every iteration later
    H = zeros(size(Z));
    T = zeros(size(Z));
    P = zeros(size(Z));
    rho = zeros(size(Z));
    a = zeros(size(Z));
    mu = zeros(size(Z));
        
% Main loop for all calculations
    for j = 1:nInput

        % Integral wants scalars so we do it separately for each altitude value
        H(j) = 1/g0 .* integral(g, 0, Z(j));

        b = 1;
        % Now that we are inside the loop, for each one of the inputs, we
        % need to identify the right layer and apply the right equations.
        % Tip: use "if" clauses, like "if H(j) <= Hb(x)". Type "help if" 
        % in the Command Window.
        if H(j) <= Hb(1)        % below sea level
            T(j) = NaN;
            P(j) = NaN;
        elseif H(j) <= Hb(2)    % b=0
            b = 1;
        elseif H(j) <= Hb(3)    % b=1
            b = 2;
        elseif H(j) <= Hb(4)    % b=2
            b = 3;
        elseif H(j) <= Hb(5)    % b=3
            b = 4;
        elseif H(j) <= Hb(6)    % b=4
            b = 5;
        elseif H(j) <= Hb(7)    % b=5
            b = 6;
        elseif H(j) <= Hb(8)    % b=6
            b = 7;
        else                    % b=7
            T(j) = NaN;
            P(j) = NaN;
        end

        b;


        T(j) = Tb(b) + bet(b) .* (H(j) - Hb(b));
        P(j) = pb(b) .* (T(j)./(T(j) + bet(b) .* (H(j) - Hb(b)))).^((g0 .* M0)/(R .* bet(b)));

        % Calculate secondary properties: density, speed of sound, and dynamic viscosity
        rho(j) = (P(j).*M0)./(T(j).*R);
        a(j) = sqrt((gamma.*R.*T(j))./(M0));
        mu(j) = (beta.*(T(j).^(1.5)))./(T(j)+S);
    
    end

% Output calculated values
%     T
%     P
%     rho
%     a
%     mu

end
