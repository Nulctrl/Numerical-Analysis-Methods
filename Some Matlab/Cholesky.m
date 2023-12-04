n=2000;
% Create a symmetric positive-definite matrix, A
A = rand(n,n);
A = A*A' + eye(n);  % Ensure A is positive definite

time_custom=0;

for i=1:5
    tic;
    C_builtin = chol(A);
    time_builtin = time_builtin+ toc;
end
time_builtin=time_builtin/5;

fprintf('Size of the System : N=%.6f \n', n);
fprintf('Time taken by custom algo: %.6f seconds\n', time_custom);
fprintf('Time taken by built in: %.6f seconds\n', time_builtin);
% Compare the results
diff = norm(C_custom - C_builtin);
disp(['Difference between custom and built-in: ', num2str(diff)]);

tolerance = 1e-10;
if diff < tolerance
    disp('The custom function is verified!');
else
    disp('The custom function did not pass the verification.');
end

function C = mychol(A)
    % Check if input matrix is square
    [n, m] = size(A);
    assert(n == m);
    C = A;  
 
    for i = 1:n
        % Check for non-positive diagonal element
        if C(i, i) <= 0
            error('Non-positive diagonal element encountered.');
        end
        % Cholesky Decomposition
        C(i, i) = sqrt(C(i, i));
        C(i, i+1:end) = C(i, i+1:end) / C(i, i);
        C(i+1:end, i) = 0;
        for j = i+1:n
            C(j, j:end) = C(j, j:end) - C(i, j:end) * C(i, j);
        end
    end
end

