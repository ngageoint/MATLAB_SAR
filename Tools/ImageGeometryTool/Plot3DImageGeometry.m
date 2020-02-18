function Vectors = Plot3DImageGeometry(meta,PlotOptions,plot_h)
%PLOTIMAGEGEOMETRY Plots Image Geometry based on meta and Plot Options

CamTarget = get(plot_h,'CameraTarget');
if CamTarget(1) == 0.5 && CamTarget(2) == 0.5 && CamTarget(3) == 0.5
    FirstPlot = 1;
else
    FirstPlot = 0;
end

Vectors = [];

%% Some Initial Computations
SCP = [meta.GeoData.SCP.ECF.X meta.GeoData.SCP.ECF.Y meta.GeoData.SCP.ECF.Z];
ARP = [meta.SCPCOA.ARPPos.X meta.SCPCOA.ARPPos.Y meta.SCPCOA.ARPPos.Z];
ARP_ = ARP./norm(ARP);
ARV = [meta.SCPCOA.ARPVel.X meta.SCPCOA.ARPVel.Y meta.SCPCOA.ARPVel.Z];
V_ = ARV./norm(ARV);

UVectECF_Col = [meta.Grid.Col.UVectECF.X meta.Grid.Col.UVectECF.Y meta.Grid.Col.UVectECF.Z];
UVectECF_Row = [meta.Grid.Row.UVectECF.X meta.Grid.Row.UVectECF.Y meta.Grid.Row.UVectECF.Z];

%SCP Geometry Structure
g = vect2geom(SCP',ARP',ARV');

R = ARP-SCP; %range vector from SCP to ARP
R_ = R./norm(R);
if isfield(meta,'PFA')
    SPN = [meta.PFA.IPN.X meta.PFA.IPN.Y meta.PFA.IPN.Z]; 
else
    SPN = cross(V_,R_);
    if meta.SCPCOA.SideOfTrack == 'L'
        SPN = -1*SPN;
    end
end

%find the SSP that is in the SCP normal plane which is the intersection of 
%the ARP vector with the SCP normal plane.  This will be at some
%elevation over the ellipsoid due to the earth curvature.
GPN = wgs_84_norm(SCP(1),SCP(2),SCP(3));
Ang = acosd(dot(GPN,ARP_));
Dist = norm(SCP)./cosd(Ang);
SSP = (ARP_)*Dist;

SSPTrue = (ARP_)*norm(SCP);

%Intersection of Ground and Slant plane
GSLine = cross(GPN,SPN);
GSLine = GSLine./norm(GSLine);

%Layover in Ground Plane
Layover = cross(GPN,GSLine);
Layover = Layover./norm(Layover);
Layover = Layover - dot(Layover,GPN').*GPN';
Layover = Layover./norm(Layover);

%Layover in Slant Plane
LayoverS = cross(GPN,GSLine);
LayoverS = LayoverS./norm(LayoverS);
LayoverS = LayoverS - dot(Layover,SPN).*SPN;
LayoverS = LayoverS./norm(LayoverS);

if meta.SCPCOA.SideOfTrack == 'L'
    CrossTrack = -1*cross(V_,GPN);
else
    CrossTrack = cross(V_,GPN);
end
CrossTrack = CrossTrack./norm(CrossTrack);
CrossTrack = CrossTrack - dot(CrossTrack,GPN').*GPN';
CrossTrack = CrossTrack./norm(CrossTrack);


%compute arc that follows earth surface between SCP and true SSP (we'll
%ignore the ellipsoid and just use a spherical earth)
NumArcPoints = 25;
NormVec = cross(GPN,ARP_);
NormVec = NormVec./norm(NormVec);
EarthArc = zeros(NumArcPoints,3);
for i=1:NumArcPoints
    RotAng = (i-1)*Ang/(NumArcPoints-1);
    T = RotateAboutAxis(NormVec,RotAng);
    EarthArc(i,:) = (T*SCP')';
end

%% Set up SCP Centered Axis and Rotate Vectors from ECF

%GPN is the Z-Axis
Z_Axis = wgs_84_norm(SCP(1),SCP(2),SCP(3))';
%project velocity into ground plane for x-axisv
X_Axis = V_ - dot(V_,Z_Axis).*Z_Axis;
%cross to Z and X to get y-axis
Y_Axis = -1*cross(X_Axis,Z_Axis);

%create ECF to SCP oriented rotation matrix
ROT = [X_Axis;Y_Axis;Z_Axis];

%Rotate all the ECF Vectors to SCP Axes.  Subtract SCP vector from anything
%other than unit vectors
ARP = ROT*(ARP'-SCP');
SSP = ROT*(SSP'-SCP');
SSPTrue = ROT*(SSPTrue'-SCP');
GPN = ROT*GPN;
SPN = ROT*SPN';
UVectECF_Col = ROT*UVectECF_Col';
UVectECF_Row = ROT*UVectECF_Row';
V_ = ROT*V_';
R_ = ROT*R_';
CrossTrack = ROT*CrossTrack';
GSLine = ROT*GSLine';
Layover = ROT*Layover';
LayoverS = ROT*LayoverS';
for i=1:NumArcPoints
    EarthArc(i,:) = (ROT*(EarthArc(i,:)'-SCP'))';
end

Vectors.SPN = SPN;
Vectors.GPN = GPN;
Vectors.R = R_;
Vectors.V = V_;
Vectors.Layover = LayoverS;
Vectors.ROT = ROT;

%Use the Range magnitude to scale the "image" and the vectors so we can see
%them at earth scales
VLength = 3*norm(R)/8;

%Ground Plane
GP(:,1) = [-1*SSP(2) SSP(2) SSP(2) -1*SSP(2) -1*SSP(2)]';
GP(:,2) = [-1*SSP(2) -1*SSP(2) SSP(2) SSP(2) -1*SSP(2)]';
GP(:,3) = [0 0 0 0 0]';

%Slant Plane
ImageRatio = double(meta.ImageData.NumCols)/double(meta.ImageData.NumRows);
SP(1,:) = -1*VLength*UVectECF_Col+VLength*UVectECF_Row/ImageRatio;
SP(2,:) = VLength*UVectECF_Col+VLength*UVectECF_Row/ImageRatio;
SP(3,:) = VLength*UVectECF_Col+-1*VLength*UVectECF_Row/ImageRatio;
SP(4,:) = -1*VLength*UVectECF_Col+-1*VLength*UVectECF_Row/ImageRatio;
SP(5,:) = -1*VLength*UVectECF_Col+VLength*UVectECF_Row/ImageRatio;

%compute Ground plane border
GPBorder = zeros(5,3);
for i=1:5
    Length = SP(i,3)/SPN(3);
    GPBorder(i,:) = SP(i,:) - SPN'*Length;
end
GPBorderLength = norm(GPBorder(2,:) - GPBorder(3,:));

%layover in slant
LayoverS = Layover - dot(Layover,SPN)*SPN;
LayoverS = LayoverS./norm(LayoverS);

%multipath in Ground Plane
Length = -1*UVectECF_Row(3)/SPN(3);
Multipath = UVectECF_Row + SPN.*Length;
Multipath = Multipath./norm(Multipath);

%Shadow in the ground plane
Shadow = -1*R_ - dot(-1*R_,GPN)*GPN;
Shadow = Shadow./norm(Shadow);

%shadow in slant plane
ShadowS = Shadow - dot(Shadow,SPN)*SPN;
ShadowS = ShadowS./norm(ShadowS);

%North in Ground
ROTZ = [cos(-1*g.azimuth) sin(-1*g.azimuth) 0; ...
        -1*sin(-1*g.azimuth) cos(-1*g.azimuth) 0; ...
        0 0 1];
North = ROTZ*SSP;    

%plot north in the slant plane
NorthS = North - dot(North,SPN)*SPN;
NorthS = NorthS./norm(NorthS);

%solve for point on the GSIntersection line that would be at broadside
VNormal = cross(SPN,V_);  %this would be the range vector at broadside
%solve for xyz point of intersection between range vector and GS
%intersection vector.  This is just a 3x3 linear system.
A = [VNormal(1) -1*GSLine(1); VNormal(2) -1*GSLine(2); VNormal(3) -1*GSLine(3)];
B = [-1*ARP(1);-1*ARP(2);-1*ARP(3)];
x = A\B;
SquintPoint = ARP + VNormal*x(1); 

%% Plot vectors and planes
hold(plot_h,'on');
LegendStrings = cell(1);

axes(plot_h);


%% set limits tight around objects
MinX = min(GP(:,1));
MaxX = max(GP(:,1));
MinY = min(GP(:,2));
MaxY = max(GP(:,2));
MinZ = min(GP(:,3));
MaxZ = max(GP(:,3));
MinX = min([MinX [ARP(1)-V_(1)*VLength/8 ARP(1)+V_(1)*VLength]]);
MaxX = max([MaxX [ARP(1)-V_(1)*VLength/8 ARP(1)+V_(1)*VLength]]);
MinY = min([MinY [ARP(2)-V_(2)*VLength/8 ARP(2)+V_(2)*VLength]]);
MaxY = max([MaxY [ARP(2)-V_(2)*VLength/8 ARP(2)+V_(2)*VLength]]);
MinZ = min([MinZ [ARP(3)-V_(3)*VLength/8 ARP(3)+V_(3)*VLength]]);
MaxZ = max([MaxZ [ARP(3)-V_(3)*VLength/8 ARP(3)+V_(3)*VLength]]);
MinX = min([MinX EarthArc(:,1)']);
MaxX = max([MaxX EarthArc(:,1)']);
MinY = min([MinY EarthArc(:,2)']);
MaxY = max([MaxY EarthArc(:,2)']);
MinZ = min([MinZ EarthArc(:,3)']);
MaxZ = max([MaxZ EarthArc(:,3)']);
MinX = min([MinX SP(:,1)']);
MaxX = max([MaxX SP(:,1)']);
MinY = min([MinY SP(:,2)']);
MaxY = max([MaxY SP(:,2)']);
MinZ = min([MinZ SP(:,3)']);
MaxZ = max([MaxZ SP(:,3)']);

%add 5% padding
MinX = MinX - round((MaxX-MinX)/20);
MaxX = MaxX + round((MaxX-MinX)/20);
MinY = MinY - round((MaxY-MinY)/20);
MaxY = MaxY + round((MaxY-MinY)/20);
MinZ = MinZ - round((MaxZ-MinZ)/20);
MaxZ = MaxZ + round((MaxZ-MinZ)/20);

set(gca,'XLim',[MinX MaxX]);
set(gca,'YLim',[MinY MaxY]);
set(gca,'ZLim',[MinZ MaxZ]);

daspect([1 1 1]);
pbaspect([1 1 1]);


%plot range line
if PlotOptions.RangeCheck
    plot3(plot_h,[0 ARP(1)],[0 ARP(2)],[0 ARP(3)],'r','LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Range Vector');
end

if PlotOptions.VelocityCheck
    plot3(plot_h,[ARP(1)-V_(1)*VLength/8 ARP(1)+V_(1)*VLength],[ARP(2)-V_(2)*VLength/8 ARP(2)+V_(2)*VLength], ...
      [ARP(3)-V_(3)*VLength/8 ARP(3)+V_(3)*VLength],'k','LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Velocity Vector');
end

if PlotOptions.GroundPlaneCheck
    patch(GP(:,1),GP(:,2),GP(:,3),'g','facecolor','g','facealpha',0.25,'edgecolor','g');
    LegendStrings = horzcat(LegendStrings,'Ground Plane');
end

if PlotOptions.SlantPlaneCheck
    patch(SP(:,1),SP(:,2),SP(:,3),'b','facecolor','b','facealpha',0.25,'edgecolor','b');
    LegendStrings = horzcat(LegendStrings,'Slant Plane');
end

if PlotOptions.GroundNormalCheck
    plot3(plot_h,[0 GPN(1)*VLength],[0 GPN(2)*VLength],[0 GPN(3)*VLength],'g','LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Ground Plane Normal');
end

if PlotOptions.SlantNormalCheck
    plot3(plot_h,[0 SPN(1)*VLength],[0 SPN(2)*VLength],[0 SPN(3)*VLength],'b','LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Slant Plane Normal');
end

if PlotOptions.GroundSlantCheck
    plot3(plot_h,[-1*GSLine(1)*VLength GSLine(1)*VLength],[-1*GSLine(2)*VLength GSLine(2)*VLength], ...
          [-1*GSLine(3)*VLength GSLine(3)*VLength],'c','LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Ground/Slant Int');
end

%need items to show up in legend if either ground or slant or both 
if PlotOptions.GroundLayoverCheck
    plot3(plot_h,[0 Layover(1)*VLength],[0 Layover(2)*VLength],[0 Layover(3)*VLength],'Color',[1 .6666 0],'LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Layover');
elseif PlotOptions.SlantLayoverCheck
    plot3(plot_h,[0 LayoverS(1)*VLength],[0 LayoverS(2)*VLength],[0 LayoverS(3)*VLength],'Color',[1 .6666 0],'LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Layover');
end

if PlotOptions.GroundMultipathCheck
    plot3(plot_h,[0 Multipath(1)*VLength],[0 Multipath(2)*VLength], ...
          [0 Multipath(3)*VLength],'Color','m','LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Multipath');
elseif PlotOptions.SlantMultipathCheck
    plot3(plot_h,[0 UVectECF_Row(1)*VLength],[0 UVectECF_Row(2)*VLength], ...
          [0 UVectECF_Row(3)*VLength],'Color','m','LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'Multipath');  
end

if PlotOptions.GroundShadowCheck
    plot3(plot_h,[0 Shadow(1)*VLength],[0 Shadow(2)*VLength], ...
          [0 Shadow(3)*VLength],'Color','k','LineWidth',2,'LineStyle','--');
    LegendStrings = horzcat(LegendStrings,'Shadow');
elseif PlotOptions.SlantShadowCheck
    plot3(plot_h,[0 ShadowS(1)*VLength],[0 ShadowS(2)*VLength], ...
          [0 ShadowS(3)*VLength],'Color','k','LineWidth',2,'LineStyle','--');
    LegendStrings = horzcat(LegendStrings,'Shadow');
end

if PlotOptions.GroundNorthCheck
    plot3(plot_h,[0 North(1)*VLength],[0 North(2)*VLength],[0 North(3)*VLength], ...
          'Color',[1 1 0],'LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'North');  
elseif PlotOptions.SlantNorthCheck
    plot3(plot_h,[0 NorthS(1)*VLength],[0 NorthS(2)*VLength],[0 NorthS(3)*VLength], ...
          'Color',[1 1 0],'LineWidth',2);
    LegendStrings = horzcat(LegendStrings,'North');  
end

if PlotOptions.GrazeCheck
    %create an arc that is 1/4 of the length of the range vector
    [Arc, ang] = ComputeArc(ARP,SSP,0.25*norm(ARP),15);  
    %fill in angle
    patch([0 Arc(:,1)'],[0 Arc(:,2)'],[0 Arc(:,3)'],'m','facecolor','m','facealpha',0.5,'edgecolor','m'); 
    a = sprintf('Graze Angle: %3.1f deg',ang);
    LegendStrings = horzcat(LegendStrings,a);  
end

if PlotOptions.AzimuthCheck
    %create an arc that is 1/4 of the length of the range vector
    if meta.SCPCOA.AzimAng > 180
        DirectionFlag = 1;
    else
        DirectionFlag = 0;
    end
    [Arc, ang] = ComputeArc(SSP,North,0.25*norm(ARP),50,DirectionFlag);  
    %fill in angle (orange)
    patch([0 Arc(:,1)'],[0 Arc(:,2)'],[0 Arc(:,3)'],'m','facecolor',[1 .6471 0],'facealpha',0.75,'edgecolor',[1 .6471 0]); 
    a = sprintf('Azimuth Angle: %4.1f deg',ang);
    LegendStrings = horzcat(LegendStrings,a);  
end
if PlotOptions.SlopeCheck
    %create an arc that is 1/4 of the length of the range vector
    [Arc, ang] = ComputeArc(GPN,SPN,0.25*norm(ARP),15);  
    %fill in angle (turquoise)
    patch([0 Arc(:,1)'],[0 Arc(:,2)'],[0 Arc(:,3)'],'m','facecolor',[0.2510 0.8784 0.8157],'facealpha',0.5,'edgecolor',[0.2510 0.8784 0.8157]); 
    a = sprintf('Slope Angle: %3.1f deg',ang);
    LegendStrings = horzcat(LegendStrings,a);       
end

if PlotOptions.TwistCheck            
    xhat = cross(SPN,ARP);
    xhat = xhat./norm(xhat);
    
    xbar = cross(GPN,SSP);
    xbar = xbar./norm(xbar);
              
    %create an arc that is 1/4 of the length of the range vector
    [Arc, ang] = ComputeArc(xbar,xhat,0.25*norm(ARP),15);  
    %fill in angle (firebrick)
    patch([0 Arc(:,1)'],[0 Arc(:,2)'],[0 Arc(:,3)'],'m','facecolor',[0.6980 0.1333 0.1333],'facealpha',0.5,'edgecolor',[0.6980 0.1333 0.1333]); 
    if meta.SCPCOA.DopplerConeAng > 90 && meta.SCPCOA.SideOfTrack == 'L'
        ang = -1*ang;
    end
    if meta.SCPCOA.DopplerConeAng < 90 && meta.SCPCOA.SideOfTrack == 'R'
        ang = -1*ang;
    end
    a = sprintf('Twist Angle: %3.1f deg',ang);
    LegendStrings = horzcat(LegendStrings,a);               
end

if PlotOptions.SquintGCheck
   %vector from SSP to origin
   A = -1*SSP;
   %vector from SSP to Squint Intersection
   B = SquintPoint - SSP;
   B = B./norm(B);
   [Arc, ang] = ComputeArc(A,B,0.25*norm(ARP),15);
   %fill in angle (orchid)
   patch([SSP(1) Arc(:,1)'+SSP(1)],[SSP(2) Arc(:,2)'+SSP(2)],[SSP(3) Arc(:,3)'+SSP(3)],'m','facecolor',[0.8549 0.4392 0.8392],'facealpha',0.5,'edgecolor',[0.8549 0.4392 0.8392]); 
   if meta.SCPCOA.DopplerConeAng > 90 %backward looking is negative
       ang = -1*ang;
   end
   a = sprintf('Ground Squint Angle: %3.1f deg',ang);
   LegendStrings = horzcat(LegendStrings,a);       
end

if PlotOptions.SquintSCheck
   %vector from SSP to origin
   A = -1*ARP;
   %vector from SSP to Squint Intersection
   B = SquintPoint - ARP;
   B = B./norm(B);
   [Arc, ang] = ComputeArc(A,B,0.25*norm(ARP),15);
   %fill in angle (dark golden rod)
   patch([ARP(1) Arc(:,1)'+ARP(1)],[ARP(2) Arc(:,2)'+ARP(2)],[ARP(3) Arc(:,3)'+ARP(3)],'m','facecolor',[0.7216 0.5255 0.0431],'facealpha',0.5,'edgecolor',[0.7216 0.5255 0.0431]); 
   if meta.SCPCOA.DopplerConeAng > 90 %backward looking is negative
       ang = -1*ang;
   end
   a = sprintf('Slant Squint Angle: %3.1f deg',ang);
   LegendStrings = horzcat(LegendStrings,a);             
end

if PlotOptions.DCACheck
    [Arc, ang] = ComputeArc(V_,-1*R_,0.25*norm(ARP),15);
    %fill in angle (dodger blue)
    patch([ARP(1) Arc(:,1)'+ARP(1)],[ARP(2) Arc(:,2)'+ARP(2)],[ARP(3) Arc(:,3)'+ARP(3)],'m','facecolor',[0.1176 0.5647 1.0000],'facealpha',0.5,'edgecolor',[0.1176 0.5647 1.0000]); 
    a = sprintf('Doppler Cone Angle: %3.1f deg',ang);
    LegendStrings = horzcat(LegendStrings,a);             
end

LegendStrings2 = [];
for i=2:length(LegendStrings)
    LegendStrings2{i-1} = LegendStrings{i};
end

%now plot objects that won't appear in legend

%SAR Point
plot3(ARP(1),ARP(2),ARP(3),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7)

%GRP Point (0,0,0)
plot3(0,0,0,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',7)

if PlotOptions.VelocityCheck
    VEnd = [ARP(1)+V_(1)*VLength ARP(2)+V_(2)*VLength ARP(3)+V_(3)*VLength];
    arrow3(ARP',VEnd);
    
    %plot cross-track vector and local horizontal
    plot3(plot_h,[ARP(1) ARP(1)+CrossTrack(1)*VLength],[ARP(2) ARP(2)+CrossTrack(2)*VLength], ...
          [ARP(3) ARP(3)+CrossTrack(3)*VLength],'k','LineWidth',1);
    plot3(plot_h,[ARP(1) ARP(1)+V_(1)*VLength],[ARP(2) ARP(2)+V_(2)*VLength], ...
          [ARP(3) ARP(3)],'k','LineWidth',1);
end

if PlotOptions.SlantLayoverCheck && PlotOptions.GroundLayoverCheck
    plot3(plot_h,[0 LayoverS(1)*VLength],[0 LayoverS(2)*VLength],[0 LayoverS(3)*VLength],'Color',[1 .6666 0],'LineWidth',2);
end
    
if PlotOptions.SlantMultipathCheck && PlotOptions.GroundMultipathCheck
    plot3(plot_h,[0 UVectECF_Row(1)*VLength],[0 UVectECF_Row(2)*VLength], ...
          [0 UVectECF_Row(3)*VLength],'Color','m','LineWidth',2);
end

if PlotOptions.SlantShadowCheck && PlotOptions.GroundShadowCheck
    plot3(plot_h,[0 ShadowS(1)*VLength],[0 ShadowS(2)*VLength], ...
          [0 ShadowS(3)*VLength],'Color','k','LineWidth',2,'LineStyle','--');
end

if PlotOptions.SlantNorthCheck && PlotOptions.GroundNorthCheck
    plot3(plot_h,[0 NorthS(1)*VLength],[0 NorthS(2)*VLength],[0 NorthS(3)*VLength], ...
          'Color',[1 1 0],'LineWidth',2);
end

if PlotOptions.RangeGroundCheck
    patch([0 ARP(1) SSP(1) 0],[0 ARP(2) SSP(2) 0],[0 ARP(3) SSP(3) 0],'r','facecolor','r','facealpha',0.25,'edgecolor','r');
end

if PlotOptions.GrazeCheck || PlotOptions.AzimuthCheck
    %local horizontal
    plot3(plot_h,[0 SSP(1)],[0 SSP(2)],[0 SSP(3)],'-k','LineWidth',1);    
end

if PlotOptions.SlopeCheck
    %also show slope between the ground and the Squint Point
    %vector from SSP to origin
    A = ARP - SquintPoint;
    A = A./norm(A);
    %vector from SSP to Squint Intersection
    B = SSP - SquintPoint;
    B = B./norm(B);
    Arc = ComputeArc(A,B,0.25*norm(ARP),15);
    %fill in angle (turquoise)
    patch([SquintPoint(1) Arc(:,1)'+SquintPoint(1)],[SquintPoint(2) Arc(:,2)'+SquintPoint(2)], ...
          [SquintPoint(3) Arc(:,3)'+SquintPoint(3)],'m','facecolor',[0.2510 0.8784 0.8157], ...
          'facealpha',0.5,'edgecolor',[0.2510 0.8784 0.8157]);         
end
    
if PlotOptions.SquintGCheck || PlotOptions.SquintSCheck || PlotOptions.SlopeCheck
    plot3(plot_h,[SSP(1) SquintPoint(1)],[SSP(2) SquintPoint(2)],[SSP(3) SquintPoint(3)],'Color','k','LineWidth',1);
    plot3(plot_h,[ARP(1) SquintPoint(1)],[ARP(2) SquintPoint(2)],[ARP(3) SquintPoint(3)],'Color','k','LineWidth',1);
    plot3(plot_h,[0 SSP(1)],[0 SSP(2)],[0 SSP(3)],'-k','LineWidth',1); 
    plot3(SquintPoint(1),SquintPoint(2),SquintPoint(3),'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6)
end

if PlotOptions.GroundBorderCheck
    plot3(GPBorder(:,1),GPBorder(:,2),GPBorder(:,3),'g','LineWidth',2);
    %Draw Lines to show the Ground Plane Trim due to squint
    %The direction will be the multipath direction
    %in the ground plane
    if meta.SCPCOA.SideOfTrack == 'L'
        if meta.SCPCOA.DopplerConeAng > 90 %backward squint   
            StartPoint1 = GPBorder(2,:);
            StartPoint2 = GPBorder(4,:);  
        else
            StartPoint1 = GPBorder(1,:);
            StartPoint2 = GPBorder(3,:);  
        end
    else
        if meta.SCPCOA.DopplerConeAng > 90 %backward squint   
            StartPoint1 = GPBorder(1,:);
            StartPoint2 = GPBorder(3,:);  
        else
            StartPoint1 = GPBorder(2,:);
            StartPoint2 = GPBorder(4,:);  
        end
    end          
        
    plot3([StartPoint1(1) StartPoint1(1)-Shadow(1)*GPBorderLength*cos(g.multipath)], ...
              [StartPoint1(2) StartPoint1(2)-Shadow(2)*GPBorderLength*cos(g.multipath)], ...
              [0 0],'Color',[0 .85 0],'LineStyle','--','LineWidth',2);
        
    plot3([StartPoint2(1) StartPoint2(1)+Shadow(1)*GPBorderLength*cos(g.multipath)], ...
              [StartPoint2(2) StartPoint2(2)+Shadow(2)*GPBorderLength*cos(g.multipath)], ...
              [0 0],'Color',[0 .85 0],'LineStyle','--','LineWidth',2);
    
end

if PlotOptions.SSPCheck
    plot3(plot_h,[ARP(1) SSPTrue(1)],[ARP(2) SSPTrue(2)],[ARP(3) SSPTrue(3)],'-k','LineWidth',1.5);
    plot3(plot_h,EarthArc(:,1),EarthArc(:,2),EarthArc(:,3),'-k','LineWidth',1.5);       
    plot3(SSP(1),SSP(2),SSP(3),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7)    
end

if ~isempty(LegendStrings2)    
    legend(plot_h,LegendStrings2);
end

set(plot_h,'XTickLabel',[]);
set(plot_h,'YTickLabel',[]);
set(plot_h,'ZTickLabel',[]);
set(plot_h,'XGrid','on');
set(plot_h,'YGrid','on');
set(plot_h,'ZGrid','on');
hold(plot_h,'off');

set(plot_h,'TickLength',[0 0]);

%% set camera position looking down the SPN if Initial Plot
if FirstPlot
    CamPos = get(gca,'CameraPosition');
    CamPos = SPN*norm(CamPos);
    set(plot_h,'CameraPosition',CamPos);
end

end

