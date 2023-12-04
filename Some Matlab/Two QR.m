[U,~] = qr(rand(30));
[V,~] = qr(rand(30));
S = diag(2 .^ (-1:-1:-30));
A = U * S * V;

[Qh, Rh] = householder(A);
[Qg, Rg] = modifiedGramSchmidt(A);

fprintf('Norm of Q^T*Q-I using Householder: %.10e\n', norm(Qh' * Qh - eye(size(Qh))));
fprintf('Norm of Q^T*Q-I using Modified Gram-Schmidt: %.10e\n', norm(Qg' * Qg - eye(size(Qg))));
fprintf('Norm of Q*R-A using Householder: %.10e\n', norm(Qh * Rh - A));
fprintf('Norm of Q*R-A using Modified Gram-Schmidt: %.10e\n', norm(Qg * Rg - A));

function [Q, R] = householder(A)
    [m, n] = size(A);
    
    if m < n
        error('m<n');
    end
    
    % Initialize
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

function [Q, R] = modifiedGramSchmidt(A)
    [m, n] = size(A);
    if m < n
        error('m<n');
    end
    % Initialize 
    Q = A;
    Q(:, 1) = Q(:, 1) / norm(Q(:, 1));

    % Orthonormalize the rest of the columns
    for i = 2:n
        for j = 1:i-1
            Q(:, i) = Q(:, i) - Q(:, j) * (Q(:, j)' * Q(:, i));
        end
        Q(:, i) = Q(:, i) / norm(Q(:, i));
    end

    R = Q' * A;
end


