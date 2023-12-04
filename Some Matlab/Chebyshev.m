f = @(x) 1 ./ (1 + 12 * x.^2);

x_fine = linspace(-1, 1, 10000);
f_true = f(x_fine);

Ns = [14,19, 24,29, 34,39, 44, 49]; 
max_errors = zeros(size(Ns)); 

for j = 1:length(Ns)
    N = Ns(j); % Current N
    k = 0:N;
    chebyshev_nodes = cos((2 * k + 1) * pi / (2 * (N + 1)));
    chebyshev_values = f(chebyshev_nodes);
    chebyshev_poly_coeffs = polyfit(chebyshev_nodes, chebyshev_values, N);
    chebyshev_poly_values = polyval(chebyshev_poly_coeffs, x_fine);
    errors = abs(f_true - chebyshev_poly_values);
    max_errors(j) = max(errors);
end

figure;
semilogy(Ns, max_errors, '-o');
xlabel('N (Number of Chebyshev points)');
ylabel('Maximal Interpolation Error');
title('Maximal Interpolation Error vs. Number of Chebyshev Points');
grid on;
