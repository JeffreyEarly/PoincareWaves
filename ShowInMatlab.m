file = '/Users/jearly/Desktop/PoincareWaves.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));

figure
quiver(x,y,u,v,0.3)
xlim([min(x) max(x)])
ylim([min(y) max(y)])