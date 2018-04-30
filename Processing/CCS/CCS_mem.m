function ccs_out = CCS_mem(ref,match)
%CCS Coherent Change Subtraction Function

% INPUTS:
%   ref          - required : complex registered reference data
%   match        - required : complex registered match data
%   
% OUTPUTS:
%   ccs_out      - Complex CCS Image
%
% VERSION:
%   1.0
%     - 20151013
%     - initial version (Based on Discussions with Ken
%                        Obenshain)
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%smoothing window...maybe add this as an input
M = 75;

%% spectral corrections
refphd = fftshift(fft2(ref));
matchphd = fftshift(fft2(match));

%phase correction
Phi = conv2(conj(refphd).*matchphd,ones(M),'same')./abs(conv2(conj(refphd).*matchphd,ones(M),'same'));
matchphd = conj(Phi).*matchphd;

%amplitude correction
RatioImage = fastrunmean(abs(refphd)./(abs(matchphd)+eps),[M M],'mean');
matchphd = matchphd.*RatioImage;

ref = ifft2(refphd);
match = ifft2(matchphd);

%% image corrections
%amplitude correction (this overall ratio seems to work best)
ratio = mean(abs(ref(:)))/mean(abs(match(:)));
ref = ref./ratio;

%phase correction (this had been judged to be unecessary)
%Phi = conv2(conj(ref).*match,ones(M),'same')./abs(conv2(conj(ref).*match,ones(M),'same'));
%match = conj(Phi).*match;

%% compute CCS
ccs_out = match-ref;

end

