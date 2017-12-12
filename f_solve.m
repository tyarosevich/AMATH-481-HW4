function x = f_solve(x0, f_coeff_a, f_coeff_b, f_coeff_c)
x0 = fft(x0);
a0 = x0(1:2048); b0 = x0(2049:4096); c0 = x0(4097:end);

a = a0.*f_coeff_a;
b = b0.*f_coeff_b;
c = c0.*f_coeff_c;

x =real(ifft([a; b; c]));
