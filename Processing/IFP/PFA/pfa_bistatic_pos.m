function [ bi_pos, freq_scale ] = pfa_bistatic_pos( tx_pos, rcv_pos, srp )
%PFA_BISTATIC_POS Compute equivalent monostatic parameters
%
% Authors: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

tx_range = tx_pos - srp;
rcv_range = rcv_pos - srp;
tx_mag = sqrt(sum(tx_range.^2,2));
rcv_mag = sqrt(sum(rcv_range.^2,2));
tx_unit = bsxfun(@rdivide, tx_range, tx_mag);
rcv_unit = bsxfun(@rdivide, rcv_range, rcv_mag);
bisector = tx_unit + rcv_unit;
bisector = bsxfun(@rdivide, bisector, sqrt(sum(bisector.^2,2))); % Unit vectors
% Equivalent bistatic position bisects tx and rcv range vectors and is at an equivalent two-way range
bi_pos = bsxfun(@times, bisector, (tx_mag + rcv_mag));
% Equivalent rf frequency is cos(beta/2)
bistatic_angle = acos(sum(tx_unit.*rcv_unit,2));
freq_scale = cos(bistatic_angle/2);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////