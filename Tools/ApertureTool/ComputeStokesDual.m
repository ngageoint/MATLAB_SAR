function [S,m,psi,chi,del] = ComputeStokesDual(HH,HV,SmoothSize)
%COMPUTESTOKESDUAL Computes Stokes parameters for dual pol complex data

%can be VV/VH as well
[ny,nx] = size(HH);

%Stokes is defines as:
% S(1): Span: HH^2+HV^2
% S(2): Difference: HH^2-HV^2
% S(3): Cov (real): 2*real(HH*HV)
% S(4): Cov (imag): -2*imag(HH*HV)
S = zeros(ny,nx,4);
S(:,:,1) = fastrunmean(abs(HH.*HH)+abs(HV.*HV),[SmoothSize SmoothSize],'mean');
S(:,:,2) = fastrunmean(abs(HH.*HH)-abs(HV.*HV),[SmoothSize SmoothSize],'mean');
S(:,:,3) = fastrunmean(2*real(HH.*conj(HV)),[SmoothSize SmoothSize],'mean');
S(:,:,4) = fastrunmean(-2*imag(HH.*conj(HV)),[SmoothSize SmoothSize],'mean');

%sum of Stokes 2-4 divided by span
m = sqrt(S(:,:,2).^2 + S(:,:,3).^2 + S(:,:,4).^2)./S(:,:,1);

%psi is the orientation angle between 2/3
psi = atan2d(S(:,:,3),S(:,:,2))/2;

chi = asind(S(:,:,4)./(m.*S(:,:,1)))/2;

%del is the phase angle between H & V
del = atan2d(S(:,:,4),S(:,:,3));

end

