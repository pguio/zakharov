function getcoefs

load linwaves7

for k=linspace(0,80,4),
	[m,ii]=min(abs(k-s.k));
	fprintf(1,'k=%8.2f nu_e=%8.2f Ek=%8.2g\n', s.k(ii), s.nue(ii), s.sEk(ii)*s.N);
end
fprintf('\n');
for k=linspace(0,80,4),
	[m,ii]=min(abs(k-s.k));
	fprintf(1,'k=%8.2f nu_i=%8.2f nk=%8.2g\n', s.k(ii), s.nui(ii), s.snk(ii)*s.N);
end
