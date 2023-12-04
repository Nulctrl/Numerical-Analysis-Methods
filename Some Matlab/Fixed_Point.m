clear all;
close all;
clc;

rng(10);
t = 2 * rand(20, 1);
b = 1 + 0.8 * t - 0.5 * t .^ 2 + 0.3 * exp(t) + 0.2 * randn(20, 1);
t(21) = 1;
b(21) = 8;

ep = 0.05;
x0 = zeros(4, 1);
count = 0;
A = ones(21, 4);
for i = 1:21
    A(i, 2) = t(i);
    A(i, 3) = t(i)^2;
    A(i, 4) = exp(t(i));
end

G = 1;
D = diag(1 ./ sqrt((A*x0-b).^2 + ep));
while G > 1e-6
    x0 = (A' * D * A) \ (A' * D * b);
    D = diag(1 ./ sqrt((A*x0-b).^2 + ep));
    G = norm(A' * D * A*x0 - A' * D * b);
    count = count + 1;
end

figure;
plot(t, b, 'bs-'); % blue squares with solid line
hold on;
fplot(@(x) 1 + 0.8 * x - 0.5 * x .^ 2 + 0.3 * exp(x), [0 2], 'g'); % green line for true model
fplot(@(x) x0(1) + x0(2)*x + x0(3)*x .^ 2 + x0(4)*exp(x), [0 2], 'r--'); % red dashed line for best fit
legend('Data points', 'True model', 'Best fit');
xlabel('t');
ylabel('y');
title('Linear Least Squares Approach with Fixed Point Method');
hold off;

CT = [];
EPS = [];
ep=0.1;
for i = 0:7
    count = 0;
    x0 = ones(4, 1);
    ep = ep / 2;
    EPS = [EPS, ep];
    G = 1;
    D = diag(1 ./ sqrt((A*x0-b).^2 + ep));
    while G > 1e-6
        x0 = (A' * D * A) \ (A' * D * b);
        D = diag(1 ./ sqrt((A*x0-b).^2 + ep));
        G = norm(A' * D * A*x0 - A' * D * b);
        count = count + 1;
    end
    CT = [CT, count];
end

figure;
plot(EPS, CT, 'bo-'); % blue circles with solid line
title('Iteration Count');
xlabel('\epsilon');
ylabel('Number of Iterations');