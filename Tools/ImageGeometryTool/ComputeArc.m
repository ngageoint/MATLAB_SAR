function [ArcData, Ang] = ComputeArc(V1,V2,Dist,NumPoints,DirectionFlag)
%COMPUTEARC Computes Arc between two 3D vectors at a specified distance

%Direction of angle, 0 indicates the shortest direction, one indicates to
%go the other direction which results in an angle > 180
if ~exist('DirectionFlag','var')
    DirectionFlag = 0;
end

%Make vectors unit vectors
V1 = V1./norm(V1);
V2 = V2./norm(V2);

%get angle between vectors
Ang = acosd(dot(V1,V2));

if DirectionFlag
    Ang = -1*(360-Ang); %make neg to go the other direction
end

%cross V1 and V2 to get axis or rotation
V3 = cross(V1,V2);
V3 = V3./norm(V3);

%allocate ArcDta
ArcData = zeros(NumPoints,3);

%Set up angles to rotate through
Angs = linspace(0,Ang,NumPoints);

for i=1:NumPoints
    %create rotation matrix to perform rotations
    T = RotateAboutAxis(V3,Angs(i));    
    ArcData(i,:) = (T*V1)'*Dist;    
end

%set angle back to positive
if DirectionFlag
    Ang = -1*Ang;
end


