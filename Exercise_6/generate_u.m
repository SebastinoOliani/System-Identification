function [u] = generate_u(n, t_max)
u=randi([1,2],n*t_max);
% shift signal to {-2;2}
for i = 1:n*t_max
    for j = 1:t_max
        if u(i,j) == 1
            u(i,j) = -2;
        end
    end
end

end

