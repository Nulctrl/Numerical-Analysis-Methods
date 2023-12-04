eig_values = wilkinson_eigenvalues();
disp('Eigenvalues:');
disp(eig_values);
roots = (10:-1:1)';
diff=eig_values-roots;
disp('Diff:');
disp(diff);

function eigenvalues = wilkinson_eigenvalues()
    roots = 1:10;
    coeffs = poly(roots);
    n = length(coeffs);
    
    A_p = zeros(n - 1); 
    for i = 2:n - 1
        A_p(i, i - 1) = 1;
    end
    for i=1:n-1
        A_p(n-i, end) = -coeffs(i+1);
    end
    eigenvalues = qr_algorithm(A_p, 1e-12);
end

function [eigenvalues] = qr_algorithm(A, tol)
    A1=A;
    max_iter = 5000;
    n = size(A, 1);
    eigenvalues = zeros(n, 1);
    iter = 0;

    residual = inf;
    
    while residual > tol && iter < max_iter
        [Q, R] = qr(A);
        [Q1, R1] = householder(A1);

        A = R * Q;
        A1 =  R1 * Q1;
        
        residual = norm(tril(A, -1), 'fro');

        iter = iter + 1;
    end
   
    eigenvalues = diag(A);
end

function [Q, R] = householder(A)
    [m, n] = size(A);
    
    if m < n
        error('m<n');
    end
    

    Q = eye(m);
    
    % Perform Householder reflections
    for j = 1:n-1
        x = A(j:m, j);
        alpha = sign(x(1)) * norm(x);
        v = x + alpha * eye(m-j+1, 1);
        v = v / norm(v);
        
        A(j:m, j:n) = A(j:m, j:n) - 2 * v * (v' * A(j:m, j:n));
        Q(j:m, :) = Q(j:m, :) - 2 * v * (v' * Q(j:m, :));
    end

    R = triu(A);
    Q = Q';
end


