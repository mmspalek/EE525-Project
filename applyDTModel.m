function [x, y] = applyDTModel(A, B, C, D, u, w)

    % keys off length of input
    
    x = zeros(height(A), length(u));
    y = zeros(height(C), length(u));
    
    for i = 1:1:length(u)
        if (i ~= length(u))
            x(:, i+1) = A*x(:, i) + B*u(:, i) + w(:, i);
        end
        y(:, i) = C*x(:, i) + D*u(:, i);
    end

end