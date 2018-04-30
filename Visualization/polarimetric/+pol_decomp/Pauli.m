function out = Pauli(data, sicd_meta)
%PAULI Pauli Decomposition
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if (nargin>1) && isfield(sicd_meta,'ImageFormation') && ...
        isfield(sicd_meta.ImageFormation,'TxRcvPolarizationProc')
    HH_ind=find(strcmpi('H:H',sicd_meta.ImageFormation.TxRcvPolarizationProc));
    HV_ind=find(strcmpi('H:V',sicd_meta.ImageFormation.TxRcvPolarizationProc));
    VH_ind=find(strcmpi('V:H',sicd_meta.ImageFormation.TxRcvPolarizationProc));
    VV_ind=find(strcmpi('V:V',sicd_meta.ImageFormation.TxRcvPolarizationProc));
else % Make band assumptions based on order
    switch size(data,3)
        case 2 % Co/cross
            HH_ind = 1; HV_ind = 2; VH_ind = []; VV_ind = [];
        case 3 % HH/cross/VV
            HH_ind = 1; HV_ind = 2; VH_ind = []; VV_ind = 3;
        case 4 % Full
            HH_ind = 1; HV_ind = 2; VH_ind = 3; VV_ind = 4;
    end
end

out = zeros(size(data,1),size(data,2),3);
% Co-pol,cross-pol
if xor(isscalar(HH_ind), isscalar(VV_ind)) && xor(isscalar(HV_ind), isscalar(VH_ind))
    %incoherent assignment (co is magenta, cross is green)
    out(:,:,1) = data(:,:,[HH_ind VV_ind]);
    out(:,:,2) = data(:,:,[HV_ind VH_ind]);
    out(:,:,3) = data(:,:,[HH_ind VV_ind]);
% HH/VV
elseif isscalar(HH_ind) && isscalar(VV_ind) && isempty(HV_ind) && isempty(VH_ind)
    out(:,:,1) = data(:,:,HH_ind) - data(:,:,VV_ind);
    out(:,:,2) = 0;
    out(:,:,3) = data(:,:,HH_ind) + data(:,:,VV_ind);
% HH,cross-pol,VV
elseif isscalar(HH_ind) && isscalar(VV_ind) && xor(isscalar(HV_ind), isscalar(VH_ind))
    out(:,:,1) = data(:,:,HH_ind) - data(:,:,VV_ind);
    out(:,:,2) = 2*data(:,:,[HV_ind VH_ind]);
    out(:,:,3) = data(:,:,HH_ind) + data(:,:,VV_ind);
% Full-pol (HH,HV,VH,VV)
elseif isscalar(HH_ind) && isscalar(VV_ind) && isscalar(HV_ind) && isscalar(VH_ind)
    out(:,:,1) = data(:,:,HH_ind) - data(:,:,VV_ind);
    out(:,:,2) = data(:,:,HV_ind) + data(:,:,VH_ind);
    out(:,:,3) = data(:,:,HH_ind) + data(:,:,VV_ind);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////