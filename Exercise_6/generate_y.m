function [y] = generate_y(u, v, a, b, n, t_max)
% system starts at rest
y=zeros(n*t_max,1);
y_i1=0;
u_i1=0;
for i=1:n*t_max
    y(i)=a*y_i1+b*u_i1+v(i);
    y_i1=y(i);
    u_i1=u(i);
end

end

