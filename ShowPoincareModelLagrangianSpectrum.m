file = '/Users/jearly/Desktop/PoincareWaves.nc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%

x = ncread(file, 'x-float');
y = ncread(file, 'y-float');
t = ncread(file, 'time');

stride = 1;
dt=t(2)-t(1);

xpos = ncread(file, 'x-position');
ypos = ncread(file, 'y-position');

% make everything a column vector
xpos = reshape(xpos, length(y)*length(x)/(stride*stride), length(t));
ypos = reshape(ypos, length(y)*length(x)/(stride*stride), length(t));

% xpos = squeeze(xpos(:,1,:));
% ypos = squeeze(ypos(:,1,:));

figure
plot(xpos(1,:),ypos(1,:))

cx = xpos' + sqrt(-1)*ypos';
cx = cx(:,1:100:end);
cv = vdiff( dt, cx, 1);

[psi,lambda]=sleptap(size(cv,1),3);
[f,spp,snn,spn]=mspec(dt,cv,psi);

% convert from cycles/second to cycles/day
f=f*86400;

figure('Position', [50 50 1200 600])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

spp(1,:) = 2*spp(1,:);
snn(1,:) = 2*snn(1,:);

plot(f,spp),ylog
hold on
plot(-f,snn)
plot(f,vmedian(spp,2), 'black', 'LineWidth', 2)
plot(-f,vmedian(snn,2), 'black', 'LineWidth', 2)

xlim([-max(f) max(f)])
vlines(-corfreq(33)*24, 'k--')