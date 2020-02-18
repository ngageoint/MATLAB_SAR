function T = RotateAboutAxis(r,theta)
%ROTATEABOUTAXIS Create a Rotation MAtrix that rotates theta through the
%unit vector r

T = zeros(3,3);

C = cosd(theta);
S = sind(theta);
t = 1-C;

T(1,1) = t*r(1)*r(1)+C;
T(1,2) = t*r(1)*r(2)-S*r(3);
T(1,3) = t*r(1)*r(3)+S*r(2);
T(2,1) = t*r(1)*r(2)+S*r(3);
T(2,2) = t*r(2)*r(2)+C;
T(2,3) = t*r(2)*r(3)-S*r(1);
T(3,1) = t*r(1)*r(3)-S*r(2);
T(3,2) = t*r(2)*r(3)+S*r(1);
T(3,3) = t*r(3)*r(3)+C;

end

