file = '/Users/jearly/Desktop/PoincareWaves.nc';
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
u3d = double(ncread(file, 'u'));
v3d = double(ncread(file, 'v'));
xpos=ncread(file,'x-position');
ypos=ncread(file,'y-position');

f0 = ncreadatt(file, '/', 'f0');
dt = t(2)-t(1)

[M, N, K] = size(u3d);

% Compute a few 'mooring' time series
cv_mooring = zeros([length(t) 1]);
subsample = 8;
iMooring = 0;
for i=1:subsample:M
	for j=1:subsample:N
		iMooring = iMooring+1;
		cv_mooring(:,iMooring) = squeeze(u3d(i,j,:) + sqrt(-1)*v3d(i,j,:));
	end
end

[fn, sn] = powspec(dt,cv_mooring);
sn = vmean(sn,2);
negativeF = find(fn<0);
positiveF = find(fn>0);
if (length(negativeF) > length(positiveF))
	negativeF(1) = [];
end
s1sided = cat(1, sn(find(fn==0)), flipud(sn(negativeF))+sn(positiveF));
f1sided = fn(find(fn>=0));
figure, plot(fn, sn), ylog
hold on, plot(f1sided, s1sided, 'k')
plot(-f1sided, s1sided, 'k')
%hold on, plot(fn, fn.^(-2), 'r')

fn = f1sided;
sn = s1sided;
% grab indices between the Coriolis and half the nyquist
fitIndices = find( fn > 1.5*f0/(2*pi) & fn < 0.5*1/(2*dt));
[P,S] = polyfit(log10(fn(fitIndices)), log10(sn(fitIndices)),1);
fit = 10^(P(2))*fn.^(P(1));
hold on
plot( fn(fitIndices), fit(fitIndices), 'r')
slope = P(1)

[psi,lambda]=sleptap(size(cv_mooring,1),3);
[fn,spp,snn,spn]=mspec(dt,cv_mooring,psi);
fn=fn/(2*pi);
figure, plot(fn, vmean(spp,2)), ylog
hold on, plot(fn, vmean(snn,2))
plot(fn,vmean(spp+snn,2))

sn = vmean(spp+snn,2);
% grab indices between the Coriolis and half the nyquist
fitIndices = find( fn > 1.5*f0/(2*pi) & fn < 0.5*1/(2*dt));
[P,S] = polyfit(log10(fn(fitIndices)), log10(sn(fitIndices)),1);
fit = 10^(P(2))*fn.^(P(1));
hold on
plot( fn(fitIndices), fit(fitIndices), 'r')
slope = P(1)


xindices = 1:2:length(x);
yindices = 1:2:length(y);
figure
iTime = 2;
quiver(x(xindices),y(yindices),u3d(xindices,yindices,iTime),v3d(xindices,yindices,iTime),0.8)
xlim([min(x) max(x)])
ylim([min(y) max(y)])