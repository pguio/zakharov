function ode1

nu = 1;
k = 2;
h = 0.01;
E0 = 1;

n = fix(5/h);
t = zeros(1,n);
En = zeros(1,n);
Ee = zeros(1,n);

En(1) = E0;
Ee(1) = E0;
for j=2:n,
  t(j) = t(j-1)+h;
  En(j) = En(j-1) - (nu+i*k^2)*En(j-1)*h;
	Ee(j) = Ee(j-1)*exp(- (nu+i*k^2)*h);
end

subplot(311),
plot(t, real(En), t, real(Ee))

subplot(312),
plot(t, imag(En), t, imag(Ee))

subplot(313)
plot(t,abs(real(En-Ee)),t,abs(imag(En-Ee)))
