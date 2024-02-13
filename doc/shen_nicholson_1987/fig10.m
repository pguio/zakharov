function fig2
% function fig2

%
% $Id: fig10.m,v 1.3 2007/05/27 14:42:49 patrick Exp $
%
% Copyright (c) 2001 Patrick Guio <patrick.guio@fys.uio.no>
%
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%

L = 2*pi;
E0 = 10.0;

k = linspace(0,20,200);

nse=growth_rate_nse(k,E0);
zak=growth_rate_zak(k,E0);

options=optimset;
knse=fminbnd('growth_rate_nse',0.0,20,options,E0,-1);
kzak=fminbnd('growth_rate_zak',0.0,20,options,E0,-1);

subplot(111)
plot(k,nse,'-',k,zak,'--');
line([knse knse],growth_rate_nse([0 knse],E0));
line([kzak kzak],growth_rate_zak([0 kzak],E0));
xlabel('wave number k')
title('growth rate \gamma(k)')
legend('nse','zakharov');

fprintf(1,'km(nse) = %f\n', knse);
fprintf(1,'km(zak) = %f\n', kzak);

orient landscape; set(gcf,'PaperOrientation','portrait');
if exist('exportfig'),
  exportfig(gcf,'fig10');
else
  print -deps fig10
end

