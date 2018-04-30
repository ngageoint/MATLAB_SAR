function ccs_out = FilterRemapCCS(CCS,RunSettings,varargin)
%FILTERREMAPCCS Applied Filtering and Remap to Complex CCS

if isempty(varargin)
    ccsmean = [];
elseif isempty(varargin{1})
    ccsmean = [];
else
    ccsmean = varargin{1};
end

ccs_out = abs(CCS);

if RunSettings.SmoothCCS
    ccs_out = fastrunmean(ccs_out,[RunSettings.CCSSmooth RunSettings.CCSSmooth],'mean');
end

if isempty(ccsmean) 
    %apply percentile cutoffs
    StartThresh = prctile(ccs_out(:),.15);
    StopThresh = prctile(ccs_out(:),99);

    ccs_out(ccs_out<StartThresh) = StartThresh;
    ccs_out(ccs_out>StopThresh) = StopThresh;

    ccs_out = ccs_out-StartThresh;
    ccs_out = ccs_out./(StopThresh-StartThresh);
else
    %mean value was passed in so just do a remap based on that...this is
    %for block processed data
    ccs_out = double(amplitudetodensity(CCS,60,40,ccsmean))./255;
end

if RunSettings.InvertCCS
    ccs_out = -1*ccs_out+1;
end

end

