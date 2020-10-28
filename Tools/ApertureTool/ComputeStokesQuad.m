function [S,m,psi] = ComputeStokesQuad(HH,HV,VH,VV,SmoothSize)
%COMPUTESTOKESQUAD: Alpha/Entropy approximation for quad pol data (used
%Stokes for dual so we're keeping the name for consistency)

%Pauli basis
P_11 = fastrunmean(abs(HH+VV).*abs(HH+VV),[SmoothSize SmoothSize],'mean');
P_22 = fastrunmean(abs(HH-VV).*abs(HH-VV),[SmoothSize SmoothSize],'mean');
P_33 = fastrunmean(abs(HV-VH).*abs(HV-VH),[SmoothSize SmoothSize],'mean');

%span is the trace of the Pauli basis
S = P_11 + P_22 + P_33;

N_11 = P_11./S;
N_22 = P_22./S;
N_33 = P_33./S;

Frob = sqrt(N_11.*N_11 + N_22.*N_22 + N_33.*N_33);
m = 1.5*(1-Frob.*Frob);
psi = acosd(sqrt(N_11));
%psi = (pi/2)*(1-N_11);


end

