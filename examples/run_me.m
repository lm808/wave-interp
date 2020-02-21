% An example file showing the interpolation of pre-computed wave data based
% upon Sharma & Dean (1981)'s second-order wave theory.
% ------------------------------------------------------------------------
% lm808, 02/2020.
% github.com/lm808, all rights reserved.

% For detailed description of the data structure, visit:
% github.com/lm808/wave-interp

clear; clc; close All

data_file = 'wave08.mat';

% extract the time vector
load(data_file, 't');

% interploate the free surface as a time-history
eta_t = fInterpEta(0, 0, t, data_file);
figure
plot(t, eta_t)
xlabel('Time [m]')
ylabel('Elevation [m]')
title('Free surface at x=0, y=0')

% interpolate the free surface as a spatial profile
xq = linspace(-65, 65, 1000);
eta_x = fInterpEta(xq, 0, 0, data_file);
figure
plot(xq, eta_x)
xlabel('x-locations [m]')
ylabel('Elevation [m]')
title('Free surface at y=0, t=0')

% interpolate the free surface as a scalar field
xq = linspace(-65, 65, 100);
yq = linspace(-65, 65, 101);
[xq, yq] = meshgrid(xq, yq);
eta_xy = fInterpEta(xq, yq, zeros(size(xq)), data_file);
figure
surf(xq, yq, eta_xy)
xlabel('x-locations [m]')
ylabel('y-locations [m]')
zlabel('Elevation [m]')
title('Free surface at t=0')

% interpolate the velocities under the crest
zq = linspace(0, max(eta_t), 1000);
u = fInterpVel(0, 0, zq, 0, data_file);
figure
plot(u, zq)
xlabel('Velocity [m/s]')
ylabel('Elevation [m]')
title('Ux under x=0, y=0 at t=0')
