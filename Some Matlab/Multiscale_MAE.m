%{
epsilon = 0.001;
x = linspace(0, 1, 1000); 

y_exact = exact_solution(epsilon, x);

Y0 = @(X, x) exp(-x) - exp(x).*exp(-X);
y_approx1 = Y0(x/epsilon, x);

y_approx2 = exp(-x) - exp(-x/epsilon);

figure;
plot(x, y_exact-y_approx1, '--r', 'LineWidth', 1); hold on;
plot(x, y_exact-y_approx2, ':b', 'LineWidth', 1);
legend('Y0 Difference ', 'In-class Difference');
xlabel('x');
ylabel('y');
title(['Difference of Solutions for \epsilon = ', num2str(epsilon)]);
grid on;


function y = exact_solution(epsilon, x)
    % Characteristic equation roots
    r = roots([epsilon, 1, 1]);
    
    % Using boundary conditions
    % y(0) = 0 => C1 + C2 = 0 => C2 = -C1
    % y(1) = 1/e => C1(e^{r1} - e^{r2}) = 1/e
    C1 = exp(-1)/(exp(r(1)) - exp(r(2)));
    C2 = -C1;

    y = C1 * exp(r(1)*x) + C2 * exp(r(2)*x);
end
%}

epsilon = 0.01; 
x = linspace(0, 1, 1000); 

y = (1/exp(1)) * (exp(x*sqrt(1-4*epsilon)/epsilon) - 1) ./ (exp(sqrt(1-4*epsilon)/epsilon) - 1) .* exp(((1-x)*(sqrt(1-4*epsilon) +1) ) / (2*epsilon));

figure;
plot(x, y, 'k', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title(['Plot of the function for \epsilon = ', num2str(epsilon)]);
grid on;




