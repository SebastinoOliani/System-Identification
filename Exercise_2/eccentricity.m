function [ret] = eccentricity(a, b, c, d, e, f)
    matrix = [a, b/2, d/2; b/2, c, e/2; d/2, e/2,f];
    eta = - sign(det(matrix));
    ret = sqrt((2*sqrt((a-c)^2+b^2)) / (eta*(a+c)+sqrt((a-c)^2+b^2)));
end