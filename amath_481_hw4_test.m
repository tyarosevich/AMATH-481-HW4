clc; clear all; close all;

eps = 1/(10^(-2)); mew = 1/(10^(-5)); f = 3; q = 2*10^(-4); 

Da = 1; Db = 1; Dc = .6;
xspan = linspace(0, 80, 2048);

L = 80; n = 2048;

k = (2*pi/L) * [0:n/2-1 (-n/2):-1];
k(1) = 10^(-7);
k2 = k.^2;
deltat = 10^(-4);
tspan = 0:deltat:.05; 

f_coeff_a = exp(Da*(-k2)*(deltat/2))';
f_coeff_b = exp(Db*(-k2)*(deltat/2))';
f_coeff_c = exp(Dc*(-k2)*(deltat/2))';

tol = 10^(-7) ; options = odeset('RelTol',tol,'AbsTol',tol);

a0 = exp(-(xspan-40).^2./20); 
x0 = [a0'; a0'; a0'];
x = x0;

for i = 1:length(tspan)
    %if mod(i, 20) == 0
    plot(xspan, x(1:2048))
    drawnow
    %end
    %t = tspan(i);
  
    x = f_solve(x, f_coeff_a, f_coeff_b, f_coeff_c);
    
    [dummyt, x] = ode45('rhs', [0 deltat/2 deltat], x, options, eps, mew, f, q);
    x = x(3,:)';
    
    x = f_solve(x, f_coeff_a, f_coeff_b, f_coeff_c);
    
end

