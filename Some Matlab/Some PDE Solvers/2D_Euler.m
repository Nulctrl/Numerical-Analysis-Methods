% 2D Euler Explicit Scheme

close all;
clear all;

dx = 0.01; dy = 0.01; dt = 0.001;
x = 0:dx:1; y = 0:dy:1;
lx = length(x); ly = length(y);

% Initial conditions
u = zeros(lx, ly);
ub = zeros(lx, ly);
u(50:55, 50:55) = 4;
ub(50:55, 50:55) = 4;
u = u';
ub = ub';

figure('Name', '2D Wave Equation - Euler Explicit Scheme', 'NumberTitle', 'off');
surfHandle = surf(x, y, u);
zlim([-1, 5]);
xlabel('x');
ylabel('y');
zlabel('u');
title('Numerical solution of 2-D wave equation');
view([-37, 80]);
colorbar;
drawnow;

% Main Iteration
for tStep = 1:500
    un = zeros(lx, ly); 

    for ix = 2:lx-1
        for iy = 2:ly-1
            un(ix, iy) = 2 * u(ix, iy) - ub(ix, iy) + (dt^2 / dx^2) * (u(ix + 1, iy) + u(ix - 1, iy) - 2 * u(ix, iy)) + (dt^2 / dy^2) * (u(ix, iy + 1) + u(ix, iy - 1) - 2 * u(ix, iy));
        end
    end

    set(surfHandle, 'ZData', un);
    drawnow;

    ub = u;
    u = un;
end
