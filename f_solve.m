function x = f_solve(x0, f_coeffs)
x0 = fft(x0);
x = ifft(x0.*f_coeffs);

