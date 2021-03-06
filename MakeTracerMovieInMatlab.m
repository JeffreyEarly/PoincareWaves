file = '/Users/jearly/Desktop/PoincareWaves2.nc';
FramesFolder ='/Users/jearly/Desktop/PoincareWavesFrames2';

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
xDomain = ncread(file, 'x');
yDomain = ncread(file, 'y');
x = ncread(file, 'x-float');
y = ncread(file, 'y-float');
t = ncread(file, 'time');

deltaX = xDomain(2)-xDomain(1);
minX = min(xDomain);
maxX = max(xDomain+deltaX);

deltaY = yDomain(2)-yDomain(1);
minY = min(yDomain);
maxY = max(yDomain+deltaY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	The stride indicates how many floats we will skip
%
stride = 4;
floatSize = 5;

% Read in the initial position of the floats.
% We will use this information to maintain a constant color on each float.
xposInitial = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) 1], [length(y)/stride length(x)/stride 1], [stride stride 1]));
xposInitial = reshape(xposInitial, length(y)*length(x)/(stride*stride), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
figure('Position', [50 50 1920 1080])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for iTime=1:length(t)
%for iTime=9:9

	% read in the position of the floats for the given time	
	xpos = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]));
	ypos = double(ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]));
	
	% make everything a column vector
	xpos = reshape(xpos, length(y)*length(x)/(stride*stride), 1);
	ypos = reshape(ypos, length(y)*length(x)/(stride*stride), 1);
	
	xpos = mod( xpos-minX, maxX-minX ) + minX;
	ypos = mod( ypos-minY, maxY-minY ) + minY;
	
	% default color map is only 128 shades---we need more!
	colormap(jet(1024))	
	
	% now plot the floats, colored by initial position
	% Scatter works, but is substantially slower than using mesh.
	% scatter(xpos, ypos, floatSize*floatSize, xposInitial, 'filled')	
	mesh([xpos';xpos'],[ypos';ypos'],[xposInitial';xposInitial'],'mesh','column','marker','.','MarkerSize',floatSize*floatSize), view(2)
	grid off
	
	% make the axes look better
	set( gca, 'TickDir', 'out');
	set( gca, 'Linewidth', 1.0);
	axis equal tight
	
	% get rid of the xticks because we're going to put a colorbar below with the same info.
	set( gca, 'xtick', [])
	
	
	xlim([minX maxX])
	ylim([minY maxY])
	
	% label everything
	title( sprintf('Floats advected by a Poincare Wave field, colored by initial position, hour %.1f', t(iTime)/3600), 'fontsize', 28, 'FontName', 'Helvetica' );
	ylabel( 'distance', 'FontSize', 24.0, 'FontName', 'Helvetica');
	
	% add a color bar
	cb = colorbar( 'location', 'SouthOutside' );
	set(get(cb,'xlabel'),'String', 'distance', 'FontSize', 24.0, 'FontName', 'Helvetica');
	set( gca, 'clim', [minX maxX] );
	
	% write everything out	
	output = sprintf('%s/%03d', FramesFolder,iTime-1);
	print('-depsc2', output)
end