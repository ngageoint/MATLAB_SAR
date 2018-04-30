function Distance = ComputeDistance( obj, pos )

meta = obj.Metadata{1};

%compute line distance and heading in the ground plane
StartPos = obj.axescoords2native(pos(1,:));
StopPos = obj.axescoords2native(pos(2,:));

DeltaX = (StartPos(1)-StopPos(1))*meta.Grid.Col.SS/cosd(meta.SCPCOA.TwistAng);
DeltaY = (StartPos(2)-StopPos(2))*meta.Grid.Row.SS/cosd(meta.SCPCOA.GrazeAng);

Distance = sqrt(DeltaX*DeltaX+DeltaY*DeltaY);

%alternate method...get lat/lon for each point and then compute the
%distance between the points
lla1 = point_slant_to_ground([StartPos(2), StartPos(1)]',meta);
lla2 = point_slant_to_ground([StopPos(2), StopPos(1)]',meta);

[d1km d2km]=lldistkm(lla1(1:2),lla2(1:2));

Distance = d1km*1000;

end