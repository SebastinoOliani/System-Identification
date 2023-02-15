function [f] = ARMAXobjective(x) % x = [theta; e]
    f = sqrt(x(7:end)'*x(7:end));
end

