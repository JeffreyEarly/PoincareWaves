file = '/Users/jearly/Desktop/PoincareWaves.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));

f0 = ncreadatt(file, '/', 'f0');

figure
iTime = 2;
quiver(x,y,u(:,:,iTime),v(:,:,iTime),0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])