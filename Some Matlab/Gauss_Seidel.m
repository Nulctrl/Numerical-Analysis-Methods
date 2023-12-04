jacobi_gauss_seidel_solver(100);
jacobi_gauss_seidel_solver(1000);
jacobi_gauss_seidel_solver(10000);


function jacobi_gauss_seidel_solver(n)
    x_j = zeros(n, 1);  
    x_gs = zeros(n, 1); 
    b = ones(n, 1);
    tol = 1e-8;
    max_iter = 1e5; 
    
    % Initial residual
    r0 = b - tridiagonal_multiply(x_j, n);
    norm_r0 = norm(r0);
    
    % Jacobi Method
    iter_j = 0;
    while true
        x_new = jacobi_iteration(x_j, b, n);
        r_j = b - tridiagonal_multiply(x_new, n);
        if norm(r_j) / norm_r0 < tol
            break;
        end
        x_j = x_new;
        iter_j = iter_j + 1;
        if iter_j > max_iter
            disp('Jacobi method did not converge within the maximum number of iterations.');
            break;
        end
    end
    
    % Gauss-Seidel Method
    iter_gs = 0;
    while true
        x_gs = gauss_seidel_iteration(x_gs, b, n);
        r_gs = b - tridiagonal_multiply(x_gs, n);
        if norm(r_gs) / norm_r0 < tol
            break;
        end
        iter_gs = iter_gs + 1;
        if iter_gs > max_iter
            disp('Gauss-Seidel method did not converge within the maximum number of iterations.');
            break;
        end
    end
    
    % Display results
    disp(['Jacobi method converged in ', num2str(iter_j), ' iterations for n=' num2str(n) '.']);
    disp(['Gauss-Seidel method converged in ', num2str(iter_gs), ' for n=' num2str(n) '.']);
end

function x_new = jacobi_iteration(x_old, b, n)
    x_new = zeros(n, 1);
    for i = 1:n
        sum = 0;
        if i > 1
            sum = sum + x_old(i-1);
        end
        if i < n
            sum = sum + x_old(i+1);
        end
        x_new(i) = (b(i) + sum) / 2.01;
    end
end

function x_new = gauss_seidel_iteration(x_old, b, n)
    x_new = x_old;
    for i = 1:n
        sum = 0;
        if i > 1
            sum = sum + x_new(i-1);
        end
        if i < n
            sum = sum + x_old(i+1);
        end
        x_new(i) = (b(i) + sum) / 2.01;
    end
end

function Ax = tridiagonal_multiply(x, n)
    Ax = zeros(n, 1);
    for i = 1:n
        Ax(i) = 2.01 * x(i);
        if i > 1
            Ax(i) = Ax(i) - x(i-1);
        end
        if i < n
            Ax(i) = Ax(i) - x(i+1);
        end
    end
end
