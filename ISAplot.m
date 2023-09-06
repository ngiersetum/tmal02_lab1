% ISA plot script - Include a brief description of your script here.   
%
% Description
%   <Description>
%
% See also
%   ISAfunction.m
%
% Requirements
%   <Function Name> (<Toolbox Name>)

% References
%   <Reference List like Papers>
%   ...
%
% Authors
%   <Author Name 1> <LiU ID 1>
%   <Author Name 2> <LiU ID 2>
%
% License
%   This program is part of an academic exercise for the course TMAL02,  
%   Linköping University, year 2023. The program/template is therefore   
%   free for non-commercial academic use.  
%
% Changes
%   <Date>  <Statement>
%   ...

%% The executable code starts here:
% Assuming that out ISA function already works and returns correct values,  
% now we should try to write some code to make a plot and visualize the  
% results.  


% Let´s start; for what values of altitude do we want to get the properties?  
% Tip: type "help colon" in the Command Window.  
altitudes = 0:50:86000;
    
% Now it might be time to send those values to our fantastic function, and  
% get back the properties we need. Tip: type "help function" in the Command  
% Window
[T, p, rho, a, mu, id1, id2] = ISAfunction(altitudes);

% By now we should have all the values we need. However, these are absolute
% values. If we want to plot all of them in the same axis it might look a
% bit unreadable. A much better idea would be to express these properties
% as fractions of their respective sea-level value. Doing this, we make 
% sure all of them are in the same order of magnitude.
T0 = 288.15; %[K]
p0 = 101325; % [Pa]
rho0 = 1.2250; %[kg/m³]
a0 = 340.294; %[m/s]
mu0 = 1.7894e-5; %[kg/ms]

TScaled = T ./ T0;
pScaled = p ./ p0;
rhoScaled = rho ./ rho0;
aScaled = a ./ a0;
muScaled = mu ./ mu0;

% Now we have what we need! Let's make a figure and plot on it. Hint: type
% "help plot" in the Command Window.
figure
hold on
grid on
xlim([0 1])

plot(TScaled(:), altitudes(:))
plot(pScaled(:), altitudes(:))
plot(rhoScaled(:), altitudes(:))
plot(aScaled(:), altitudes(:))
plot(muScaled(:), altitudes(:))

% And, in order to make the figure look decent, we can add and change some
% of the plot properties (subscripts are like in Latex Math mode indicated
% by "_{myText}"):
xlabel('T/T_{0} , P/P_{0} , \rho/\rho_{0} , a/a_{0} , \mu/\mu_{0} [-]');
ylabel('Z [km]');
Title('...');
legend('T/T_{0}', 'P/P_{0}', '\rho/\rho_{0}', 'a/a_{0}', '\mu/\mu_{0}', 'Location', 'NorthWest');
...

