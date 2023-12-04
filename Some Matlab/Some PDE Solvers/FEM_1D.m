function FEM_1D
    result = F_E_solver_1D_linear();
    disp(result);

    function result = FE_basis_local_1D(x, vertices, basis_type, basis_index, derivative_order)
        xnp1 = vertices(2);
        xn = vertices(1);
        h = xnp1 - xn;

        if basis_type == 1
            if basis_index == 1
                if derivative_order == 0
                    result = (xnp1 - x) / h;
                elseif derivative_order == 1
                    result = -1 / h;
                else
                    result = 0;
                end
            elseif basis_index == 2
                if derivative_order == 0
                    result = (x - xn) / h;
                elseif derivative_order == 1
                    result = 1 / h;
                else
                    result = 0;
                end
            else
                error('Invalid basis index!');
            end
        else
            error('Invalid basis type!');
        end
    end

    function result = F_E_solver_1D_linear()
        [P, T] = generate_mesh(0, 1, 0.01); 

        A = assemble_matrix_1D_linear(P, T); 
        b = assemble_vector_1D_linear(P, T); 
        [A, b] = treat_Dirichlet_boundary_1D(A, b, P, T);

        result = A \ b;
    end

    function int_value = Gauss_quatrature_1D_linear_trial_test(coe_fun, Gauss_weight, Gauss_node, vertices, basis_type_trail, basis_index_trail, derivative_order_trail)
        Gpn = length(Gauss_weight);
        int_value = 0;
        for k = 1:Gpn
            int_value = int_value + Gauss_weight(k) * feval(coe_fun, Gauss_node(k)) * FE_basis_local_1D(Gauss_node(k), vertices, basis_type_trail, basis_index_trail, derivative_order_trail) * FE_basis_local_1D(Gauss_node(k), vertices, basis_type_trail, basis_index_trail, derivative_order_trail);
        end
    end

    function A = assemble_matrix_1D_linear(P, T)
        matrix_size = [length(P), length(P)];
        number_of_elements = size(T, 2);
        Nlb = 2;
        Tb = T;  

        A = zeros(matrix_size(1), matrix_size(2));

        for n = 1:number_of_elements
            for alpha = 1:Nlb
                for beta = 1:Nlb
                    int_value = Gauss_quatrature_1D_linear_trial_test(fp);
                    A(Tb(beta, n), Tb(alpha, n)) = A(Tb(beta, n), Tb(alpha, n)) + int_value;
                end
            end
        end
    end

    function [P, T] = generate_mesh(left, right, h)
        N = (right - left) / h;
        P = left:h:right;
        T = zeros(2, N);
        for i = 1:N
            T(:, i) = [i, i + 1]';
        end
    end
end