n = 5000;
L = tril(rand(n, n),-1) + eye(n); 
b = rand(n, 1);    

% Verification and Timing
tic;
x1 = forward1(L, b);
time_forward1 = toc;

tic;
x2 = forward2(L, b);
time_forward2 = toc;

tic;
x3 = forward3(L, b);
time_forward3 = toc;

tic;
x_builtin = L \ b;
time_builtin = toc;

fprintf('Time taken by forward1: %.6f seconds\n', time_forward1);
fprintf('Time taken by forward2: %.6f seconds\n', time_forward2);
fprintf('Time taken by forward3: %.6f seconds\n', time_forward3);
fprintf('Time taken by built-in: %.6f seconds\n', time_builtin);


function x = forward1(L, b)
    % Check if conditions are met
    [n, m] = size(L);
    assert(n==m);
    assert(n==length(b));
    assert(istril(L));

    x = zeros(n, 1);
    % Forward substitution
    for i = 1:n
        sum = 0;
        for j = 1:i-1
            sum = sum + L(i, j) * x(j);
        end
        x(i) = (b(i) - sum) / L(i, i);
    end
end

function x = forward2(L, b)
    % Check if conditions are met
    [n, m] = size(L);
    assert(n==m);
    assert(n==length(b));
    assert(istril(L));

    x = zeros(n, 1);
    % Forward substitution with vector operations
    x(1) = b(1) / L(1, 1);  % Solving for the first x directly  
    for i = 2:n
        x(i) = (b(i) - L(i, 1:i-1) * x(1:i-1)) / L(i, i);
    end
end

function x = forward3(L, b)
    % Check if conditions are met
    [n, m] = size(L);
    assert(n==m);
    assert(n==length(b));
    assert(istril(L));
    
    % Forward substitution with vector operations in columns
    b(1) = b(1) / L(1, 1);
    for j = 2:n
        b(j) = b(j) - L(j, 1:j-1) * b(1:j-1);
        b(j) = b(j) / L(j, j);
    end
    x = b;
end

