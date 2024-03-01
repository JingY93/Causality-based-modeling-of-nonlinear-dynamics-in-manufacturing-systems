function [ x, y,z,a] = coupled_logistic( x1, y1,z1,a1, N )
%COUPLED_LOGISTIC(x1, y1, rx, ry, bxy, byx, N) Computes N iteration of the
% coupled logistic map, with initial values x1, y1; control parameters
% (growth rates) rx, ry and coupling between the variables bxy and byx.
x = zeros(1,N);
y = zeros(1,N);
z = zeros(1,N);
a = zeros(1,N);
x(1) = x1;
y(1) = y1;
z(1) = z1;
a(1) = a1;
for t = 2:N
    x(t) = x(t-1)*(3.9-3.9*x(t-1));
    y(t) = y(t-1)*(3.6-0.4*x(t-1)-3.6*y(t-1));
    z(t) = z(t-1)*(3.6-0.4*y(t-1)-3.6*z(t-1));
    a(t) = a(t-1)*(3.8-0.35*z(t-1)-3.8*a(t-1));
end

end