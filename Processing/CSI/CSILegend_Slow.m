function CSILegend_Slow(filename,energy_down,varargin)
%plot polar diagram corresponding to CSI
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

reader_obj = open_reader(filename);
if iscell(reader_obj)
    reader_obj = reader_obj{1};
end
meta = reader_obj.get_meta();
reader_obj.close();

try % Metadata may not be available to compute azimuth angles.
    NumLines = 512;
    AzRange = ComputeAz(meta,100*(0:NumLines)/NumLines,varargin{:});
catch % Can't create figure.
    warning('CSILegend_Slow:InsufficientMetadata','Insufficient metadata to create CSI legend.');
    return;
end

if (energy_down)
    Azimuth = ComputeAz(meta,50,varargin{:});
else % north up
    Azimuth = 0;
end

%set up polar map
set(figure,'Position',[100 100 500 500]);
hold on;
UnitCircleX = cosd(1:360);
UnitCircleY = sind(1:360);
plot(UnitCircleX,UnitCircleY,'-k','LineWidth',2);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'Box','on');
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

%plot major/minor heading lines
North = 90 - 0   + Azimuth;
South = 90 - 180 + Azimuth; 
East = 90 - 90 + Azimuth; 
West = 90 - 270 + Azimuth; 
%major
x = [cosd(North),cosd(South)];
y = [sind(North),sind(South)];
plot(x,y,'-k','LineWidth',1);
x = [cosd(East),cosd(West)];
y = [sind(East),sind(West)];
plot(x,y,'-k','LineWidth',1);
%minor
x = [cosd(North+45),cosd(South+45)];
y = [sind(North+45),sind(South+45)];
plot(x,y,'--k','LineWidth',1);
x = [cosd(East+45),cosd(West+45)];
y = [sind(East+45),sind(West+45)];
plot(x,y,'--k','LineWidth',1);

%label N,S,E,W
text(1.1*cosd(North),1.1*sind(North),'N');
text(1.1*cosd(South),1.1*sind(South),'S');
text(1.1*cosd(East),1.1*sind(East),'E');
text(1.1*cosd(West),1.1*sind(West),'W');

%for each angle plot the corresponding jet color
cmap = flipud(jet(length(AzRange)));
for i=1:length(AzRange)
    Angle = 90 - (AzRange(i) - Azimuth);
    color = cmap(i,:);
    x = [0,cosd(Angle)];
    y = [0,sind(Angle)];
    plot(x,y,'Color',color);
end

%label figure
IID = meta.CollectionInfo.CoreName(1:min(16,end));
title(sprintf('Slow-Time CSI Legend for IID: %s\n StartAz: %4.1f, StopAz: %4.1f',...
    IID, AzRange(1), AzRange(end)));

hold off;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////