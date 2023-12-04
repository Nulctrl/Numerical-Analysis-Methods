L = [0 0 1 0.5; 
     1/3 0 0 0; 
     1/3 0.5 0 0.5; 
     1/3 0.5 0 0];
v = [1 1 1 1]';
[eigenval, eigenvec] = eigen(L, v);

fprintf('Eigenvalue (eigenval) : %f\n', eigenval);
fprintf('Eigenvector (eigenvec) :\n');
disp(eigenvec);

function [a, b] = eigen(L, v)
    [m, n] = size(L);
    if m ~= n
        error('nonsquare');
    end
    
    err = inf;
    eigen = 0;
    
    while err > 1e-10
        v = L * v / norm(L * v);
        eigen = v' * L * v;
        err = norm(L * v - eigen * v) / norm(v);
    end
    
    a = eigen;
    b = v;
end
