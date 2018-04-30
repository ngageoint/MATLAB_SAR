function out = pedfremap(x)
% PEDFREMAP Piecewise extended density format remap
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

    out = amplitudetodensity(x);
    out(out > 128) = .5 * (out(out>128) + 128);
    out = uint8(out);
end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
