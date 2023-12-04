% Plotting l_k(x) for k from 0 to 7
X = linspace(-0.5, 5, 1000); % Define the range of x
hold on;
for k = 0:7
    Y = my_lk(k, X);
    plot(X, Y);
end
hold off;

% Label the axes and add a legend
xlabel('x');
ylabel('$$l_k(x)$$','Interpreter', 'latex');
title('Plot of $$l_k(x)$$ for k from 0 to 7','Interpreter', 'latex');
grid on;
legend(arrayfun(@(k) ['l_', num2str(k), '(x)'], 0:7, 'UniformOutput', false));

function lk = my_lk(k, X) 
    % initialization
    lk1 = ones(size(X)); % value of l_0(x)
    lk2 = 1 - X; % value of l_1(x)
    
    if k == 0
        lk = lk1;
    elseif k == 1
        lk = lk2;
    else
        for i = 1:k-1
            lk3 = ((2*i+1-X).*lk2/(i+1)) - i*lk1/(i+1);
            lk1 = lk2;
            lk2 = lk3;
        end
        lk = lk2;
    end
end

