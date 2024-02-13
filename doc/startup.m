%startup.m

%set(0,'DefaultFigurePosition',[550 500 505 405]);
set(0,'DefaultFigurePaperType','a4letter');
set(0,'DefaultFigurePaperUnits','centimeters');
%set(0,'DefaultFigurePaperPosition',[2.5 2.5 16.5 24.2]/2.53807);
set(0,'DefaultFigurePaperPosition',[2.5 2.5 16.5 24.2]);
set(0,'DefaultFigureNumberTitle','off');

set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontWeight','normal');
set(0,'DefaultAxesFontSize',14);
set(0,'DefaultAxesTickDir','in');
set(0,'DefaultAxesXGrid','on');
set(0,'DefaultAxesYGrid','on');
set(0,'DefaultAxesZGrid','on');

set(0,'DefaultTextFontName','times');
set(0,'DefaultTextFontWeight','normal');
set(0,'DefaultTextFontSize',16);

set(0,'DefaultLineMarkersize',2);
set(0,'DefaultLineLineWidth',.2);

set(0,'DefaultSurfaceLineWidth',.65);

addpath('../matlab')

if exist('sd')~=2,
	addpath([getenv('HOME') '/research/codes'])
end

