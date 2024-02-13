function fig48

a = linspace(1e-2,1e1,100);

TeoverTi = [0.5 1 2 3];
Fi = zeros(100,4);
for i=1:4,
	Fi(:,i) = 1./(1+a.^2)./(1+a.^2+TeoverTi(i));
end
Fe = a.^2./(1+a.^2);

semilogx(a,Fe, a, Fi);
