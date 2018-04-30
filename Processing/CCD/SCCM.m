function [sccm, anglemap] = SCCM( reference_image, match_image, corr_window_size,NR,NM )
%SCCM Signal Correlation Change Metric
%    [sccm, anglemap] = SCCM( reference_image, match_image, corr_window_size,NR,NM )
%
% Computes Signal Correlation Change Metric, a process described by
% the briefing "An Improved CCD Metric" presented by Tom Mitchell at the
% 2016 March Tech Forum
%
% Written by: Tom Mitchell, NGA-W/AIM
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[ccd, anglemap] = ccdmem(reference_image,match_image,corr_window_size);

RstarM = conv2(conj(reference_image).*match_image, ones(corr_window_size), 'same');
MstarR = conv2(conj(match_image).*reference_image, ones(corr_window_size), 'same');

R_squared = conv2(conj(reference_image).*reference_image, ones(corr_window_size), 'same');
M_squared = conv2(conj(match_image).*match_image, ones(corr_window_size), 'same');

sigR=R_squared./(2*corr_window_size^2)-NR;
sigM=M_squared./(2*corr_window_size^2)-NM;

sccm = (RstarM+MstarR)./(4*corr_window_size^2.*sqrt(sigR.*sigM));

sccm(sigR<=0)=-100; % set no signal pixels to -100
sccm(sigM<=0)=-100; % set no signal pixels to -100


end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////