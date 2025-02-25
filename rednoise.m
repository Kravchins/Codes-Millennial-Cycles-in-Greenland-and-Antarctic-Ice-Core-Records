function [noise] = rednoise(m)

x = randn(1, m);
X = fft(x); % Calculate FFT
u = m/2 + 1;
k = 1:u;
X = X(1:u);      
X = X./k;
X = [X conj(X(end-1:-1:2))];
y = real(ifft(X));
y = y(1, 1:m);
y = reshape(y, [m,1]);
noise=(y-mean(y))./std(y);

end