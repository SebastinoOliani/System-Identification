function [c,ceq] = ARMAXconstraint(x, y, PhiTyu, K)
e = x(7:end);
theta = x(1:6);
PhiTe = zeros(K,2);
PhiTe(2,1) = e(1);

for j = 3:K
    PhiTe(j,:) = [e(j-1), e(j-2)];
end

ceq = y - [PhiTyu, PhiTe] * theta - e; 
c = [];

end

