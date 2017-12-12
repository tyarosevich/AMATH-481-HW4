function x = rhs(t, x0, dummy, eps, mew, f, q)

a0 = x0(1:2048); b0 = x0(2049:4096); c0 = x0(4097:end);
ab0 = a0.*b0;

da = mew*(-q*a0 - ab0 + f*c0);
db = eps * (q*a0 - ab0 + b0.*(1 - b0));
dc = b0 - c0;

x = [da; db; dc];
