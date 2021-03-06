clc; clear all; close all;

eps = 1/(10^(-2)); mew = 1/(10^(-5)); f = 3; q = 2*10^(-4); 
L = 80; n = 2048;
Da = 1; Db = 1; Dc = 0.6;
xspan = linspace(0, 80, n)';



k = (2*pi/L) * [0:(n/2-1) (-n/2):-1];
k(1) = 10^(-7);
k2 = (k.^2)';
deltat = 10^(-4);
tspan = 0:deltat:0.05; 

f_coeff_a = exp(Da*(-k2)*(deltat/2));
f_coeff_b = exp(Db*(-k2)*(deltat/2));
f_coeff_c = exp(Dc*(-k2)*(deltat/2));
f_coeffs = [f_coeff_a; f_coeff_b; f_coeff_c];

tol = 10^(-7) ; options = odeset('RelTol',tol,'AbsTol',tol);

a0 = exp(-((xspan-40).^2)/20); 
x0 = [a0; a0; a0];
x = x0;

for j = 1:length(tspan)
    subplot(2, 2, 1)
    plot(xspan, x(1:2048))
    title('a(x)')
    subplot(2, 2, 2)
    plot(xspan, x(2049:4096))
    title('b(x)')
    subplot(2, 2, 3)
    plot(xspan, x(4097:end))
    title('c(x)')
    drawnow
    %t = tspan(i);
    
    x = f_solve(x, f_coeffs);
    x = max(0, x);
    
    [dummyt, x] = ode45('rhs', linspace(0, deltat, 3), x, options, eps, mew, f, q);
    x = x(2,:)';
    x = max(0, x);
    
    x = f_solve(x, f_coeffs);
    x = max(0, x);
    
end

%% HW 4 Part 2
clc; clear all; close all;
tic
n = 50;
Lx = 10; Ly = 12; Nx = 20; Ny = 24;
beta = 1; Diff1 = 0.1; Diff2 = Diff1;
tspan = 0:0.2:4;

[D, x] = cheb(n);
y = x;


D2 = D^2;
D2(1,:) = 0; D2(end, :) = 0;
I = eye(n + 1);
DD2 = kron(I, 100 * D2) + kron(144 * D2, I);

[X, Y] = meshgrid(x, y);

u = tanh(sqrt(X.^2+Y.^2)).*cos(angle(X+i*Y)-(sqrt(X.^2+Y.^2))) ;
v = tanh(sqrt(X.^2+Y.^2)).*sin(angle(X+i*Y)-(sqrt(X.^2+Y.^2))) ;
u = reshape(u, (n+1)^2, 1);
v = reshape(v, (n+1)^2, 1);
y0 = [u; v];
tol = 10^(-7) ; options = odeset('RelTol',tol,'AbsTol',tol);

[t, ysol] = ode45('rhs_cheb', tspan, y0, options, DD2, Diff1, Diff2);
toc

