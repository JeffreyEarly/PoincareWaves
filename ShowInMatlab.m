file = '/Users/jearly/Desktop/PoincareWaves.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));

f0 = ncreadatt(file, '/', 'f0');
T = 2*pi/f0;
inertialIndex = find(t>T, 1, 'first')

figure
subplot(2,2,1)
quiver(x,y,u(:,:,1),v(:,:,1),0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])

subplot(2,2,2)
quiver(x,y,u(:,:,floor(inertialIndex/4)),v(:,:,floor(inertialIndex/3)),0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])

subplot(2,2,3)
quiver(x,y,u(:,:,floor(2*inertialIndex/4)),v(:,:,floor(2*inertialIndex/4)),0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])

subplot(2,2,4)
quiver(x,y,u(:,:,floor(3*inertialIndex/4)),v(:,:,floor(3*inertialIndex/4)),0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])


figure
subplot(1,2,1)
plot(x, u(1,:,1), 'b')
hold on
plot(x, v(1,:,1), 'r')
subplot(1,2,2)
plot(y, u(:,1,1), 'b')
hold on
plot(y, v(:,1,1), 'r')

xpos = ncread(file, 'x-position');
ypos = ncread(file, 'y-position');

figure, plot(squeeze(xpos(1,1,:)), squeeze(ypos(1,1,:)))
dt=t(2)-t(1);
u1 = vdiff(dt, squeeze(xpos(1,1,:)),1);
v1 = vdiff(dt, squeeze(ypos(1,1,:)),1);
figure, plot(t, sqrt(u1.*u1 + v1.*v1))

figure, scatter(t, squeeze(xpos(1,1,:)))
figure, scatter(t, squeeze(ypos(1,1,:)))