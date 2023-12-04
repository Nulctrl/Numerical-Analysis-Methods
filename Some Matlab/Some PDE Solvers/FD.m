%Finite Difference Scheme
clear all;
close all;

% Initialize parameters
xl = 0; xr = 2;
M = 200;
dx = (xr-xl) / M;
tf = 2;
Nt = 800;
dt = tf/Nt;
c1 = 50;
c = 1;
mu = c^2 * dt / dx;

% Initial conditions
x = xl:dx:xr;
f = exp(-c1*(x-1).^2);
g = zeros(1, M+1); % Assuming zero initial velocity
u = zeros(Nt, M+1);
u(1, :) = f; 
udata = f;
tdata = 0;

fig = figure(1);
ax = axes('Parent', fig);
h = animatedline(ax, 'Color', 'b');
axis([xl xr -0.5 1.5]);
xlabel(ax, 'x');
ylabel(ax, 'u');
title(ax, 't=0.00');
box on;

set(fig, 'CloseRequestFcn', @(src, event)delete(fig));

% Main Iteration
for n = 1:Nt
    if ~isvalid(fig)
        break;
    end
    
    if n == 1
        u(n, 1) = u(n, M);        
        u(n, M+1) = u(n, 2);      
    elseif n == 2
        for m = 2:M
            u(n, m) = dt * g(m) + u(n - 1, m);
        end
        u(n, 1) = u(n, M);        
        u(n, M+1) = u(n, 2);      
    else
        for m = 2:M
            u(n, m) = mu^2 * (u(n - 1, m + 1) + u(n - 1, m - 1)) + 2 * (1 - mu^2) * u(n - 1, m) - u(n - 2, m);
        end
        u(n, 1) = u(n, M);        
        u(n, M+1) = u(n, 2);      
    end

    clearpoints(h);
    addpoints(h, x, u(n, :));
    title(ax, sprintf('t =%5.2f', n * dt));
    drawnow;

    if mod(n, 100) == 0
        udata = [udata; u(n, :)];
        tdata = [tdata; n * dt];
    end
end

if isvalid(fig)
    figure(2);
    waterfall(x, tdata, udata);
    view(5, 40);
    axis([xl xr 0 tf -0.5 3]);
    xlabel('x');
    ylabel('t');
    zlabel('u');
end
