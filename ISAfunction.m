function [T, P, rho, a, mu, liuID1, liuID2] = ISAfunction(Z)
% ISA function - Calculate atmospheric properties for one or several
% altitudes
%
% Synopsis
%   [T P rho a mu liuid1 liuid2] = ISAfunction(Z)
%
% Description
%   <Description>
%
% Inputs ([]s are optional)
%   (scalar or vector) Z      Geometric altitude in [m]
%
% Outputs ([]s are optional)
%   (scalar or vector) T      Temperature, K
%   (scalar or vector) P      Pressure, Pa
%   (scalar or vector) rho    Density, kg/m^3
%   (scalar or vector) a      Sound speed, m/s
%   (scalar or vector) mu     Dynamic viscosity, kg/(m s)
%
% Examples
%   <Example Code>
%
% See also
%   <See Also Function Name>
%
% Requirements
%   <Function Name> (<Toolbox Name>)

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
% Changes
%   <Date>  <Statement>
%   ...

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
    bet = [-6.5, 0, 1, 2.8, 0, -2.8, -2.0];
    Tb =[288.15, 216.650, 216.650, 228.650, 270.650, 270.560, 214.650, 186.8673]; %[K]
    pb = [101325, 22632, 5474.8, 868.01, 110.9, 66.938, 3.9564, 0.39814]; % [Pa]
    rhob = [1.2250]; %[kg/m³]
    M0 = 28.9644; %[kg/kmol]
    ab = [340.294]; %[m/s]
    mub = [1.7894e-5]; %[kg/ms]
    
    
% Determine number of input altitudes
    nInput = length(Z);
    
% Convert input (geometrical altitude) to geopotential altitude for
% calculations
    g = @(Z) g0 .* (r0./(r0+Z)).^2;   % function for calculation of g
%     H = 1/g0 .* integral(g, 0, Z);
    H = Z;   % initialize array just so it doesn't change size every iteration later
        
% What are the "primary" properties that we need to compute from the
% reference values and their equations? Which others can be computed later
% from the primary ones? We start then with the "primary" ones. If instead
% of a single value, the input contains multiple values, we will need to
% make a loop so the equations are applied to each one of those inputs.
% Let's make a loop as long as the number of inputs we caounted before!
% Type "help for" in the Command Window".
    for j = 1:nInput

        H(j) = 1/g0 .* integral(g, 0, Z(j))   % integral wants scalars so we do it separately for each altitude value

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

        b


        T(j) = Tb(b) + bet(b) .* (H(j) - Hb(b))
        P(j) = pb(b) .* (T(j)./(T(j) + bet(j) .* (H(j) - Hb(b)))).^((g0 .* M0)/(R .* bet(b)))

        % Calculate secondary properties: density, speed of sound, and dynamic viscosity
        rho(j) = (P(j).*M0)./(T(j).*R);
        a(j) = sqrt((gamma.*R.*T(j))./(M0));
        mu(j) = (beta.*(T(j).^(1.5)))./(T(j)+S);
    
    end

% Output calculated values
    T
    P
    rho
    a
    mu

end

% function g = calculate_g(Z)
%     r0 = 6356766;
%     g0 = 9.80665;
%     g = g0 + (r0/(r0+Z))^2;
% end
