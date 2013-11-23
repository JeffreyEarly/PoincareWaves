file = '/Users/jearly/Desktop/PoincareWaves.nc';
FramesFolder ='/Users/jearly/Desktop/PoincareWavesVelocityfieldFrames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(FramesFolder) == 0
	mkdir(FramesFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 'time');
u = double(ncread(file, 'u'));
v = double(ncread(file, 'v'));
xpos = ncread(file, 'x-position');
ypos = ncread(file, 'y-position');

minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);

xpos = mod( xpos-minX, maxX-minX ) + minX;
ypos = mod( ypos-minY, maxY-minY ) + minY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
figure('Position', [50 50 1920 1080])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for iTime=1:length(t)

	quiver(x,y,u(:,:,iTime),v(:,:,iTime),0.8)
	xlim([min(x) max(x)])
	ylim([min(y) max(y)])
	
	hold on
	scatter(xpos(1,1,iTime), ypos(1,1,iTime),10^2,'MarkerEdgeColor','b','MarkerFaceColor','c','LineWidth',1.5)
	hold off
	
	% write everything out	
	output = sprintf('%s/Day_%03d', FramesFolder,iTime-1);
	print('-depsc2', output)
end