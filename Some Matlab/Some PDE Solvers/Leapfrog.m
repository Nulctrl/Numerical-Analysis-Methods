% Leap-frog Scheme

clear all;
close all;

% Initialize parameters
xl = 0; xr = 2;
M = 400;
dx = (xr-xl) / M;
tf = 2.0;
Nt = 800;
dt = tf/Nt;
c1 = 50;
c = 1;
mu = c * dt / dx;

% Initial conditions
x = xl:dx:xr;
f = exp(-c1*(x-0.3).^2);
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

% Main Itration
for n = 2:Nt
    if ~isvalid(fig)
        break;
    end

    if n == 2 % 
        for m = 2:M
            u(n, m) = u(n-1, m) - mu * (u(n-1, m+1) - u(n-1, m-1));
        end
    else  
        for m = 2:M
            u(n, m) = u(n-2, m) - mu * (u(n-1, m+1) - u(n-1, m-1));
        end
    end
    u(n, 1) = u(n, M);
    u(n, M+1) = u(n, 2);

    clearpoints(h);
    addpoints(h, x, u(n, :));
    title(ax, sprintf('t =% 5.2f', n * dt));
    drawnow;

    if mod(n, 100) == 0
        udata = [udata; u(n, :)];
        tdata = [tdata; n * dt];
    end
end

if isvalid(fig)
    figure(2);
    waterfall(x, tdata, udata);
    view(10, 60);
    axis([xl xr 0 tf -1 3]);
    xlabel('x');
    ylabel('t');
    zlabel('u');
end
