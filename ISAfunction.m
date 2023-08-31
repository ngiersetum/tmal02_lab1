function [T, P, rho, a, mu, liuID1, liuID2] = ISAfunction(Z)
% ISA function - Include a brief description here.   
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
  liuID1 = "examp978" % author 1 name / <LiU ID 1>
  liuID2 = "namex123" % author 2 name / <LiU ID 2> (leaf empty [] if one person only)
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
    bs = [0, 1, 2, 3, 4, 5, 6, 7];
    S = 110; %[K]
    Hb =[0, 11000, 20000, 32000, 47000, 51000, 71000, 84852]; %[m]
    bet = [-6.5, 0, 1, 2.8, 0, -2.8, -2.0];
    Tb =[reference temperature values];
    pb = [reference pressure values];
    rhob = 1.2250; %[kg/m³]
    ab = 340.294; %[m/s]
    mub = 1.7894e-5; %[kg/ms]
    
    
% How many input values do we have? One? Two? Or several? Maybe we need to 
% count them! Tip: type "help length" in the Command Window.
    nInput = ...
    
% The input values come in geometrical altitude, right? Do we need to
% convert them? How?
    ...
        
% What are the "primary" properties that we need to compute from the
% reference values and their equations? Which others can be computed later
% from the primary ones? We start then with the "primary" ones. If instead
% of a single value, the input contains multiple values, we will need to
% make a loop so the equations are applied to each one of those inputs.
% Let's make a loop as long as the number of inputs we caounted before!
% Type "help for" in the Command Window".
    for j = 1:nInput
        
        % Now that we are inside the loop, for each one of the inputs, we
        % need to identify the right layer and apply the right equations.
        % Tip: use "if" clauses, like "if H(j) <= Hb(x)". Type "help if" 
        % in the Command Window.
        if Z(j) <= Hb(x)
                X(j) = ...   
        elseif ...
                X(j) = ...   
        else
                ...
        end
    
        % We are still inside the main loop. Now that we have calculated 
        % the "primary" properties for this step, we can use them to 
        % calculate the rest of the properties, right?
        X(j) = ...;
        X(j) = ...;
        X(j) = ...;
        ...
    
    end

% It seems we are done now! The output variables will be returned from the
% function. If you want, you can name the variables here so they will be
% displayed in the command window:
    X
    Y
    ...

end
