function dx = rhs_cheb(t, x0, dummy, DD2, Diff1, Diff2)
tic
u = x0(1:2601); v = x0(2602:end);

Lu = DD2*u;
Lv = DD2*v;
A2 = (u.^2) .* (v.^2);
lambda = 1 - A2;
omega = -A2;

du = lambda .* u - omega .* v + Diff1*Lu;
dv = omega .* u + lambda .* v + Diff2*Lv;

dx = [du; dv];
toc