function [ output_meta ] = meta2sicd_s1cal( domnode )
%META2SICD_S1CAL Converts Sentinel-1 calibration XML file into SICD format
%
% Written by: Wade Schwartzkopf, NGA Research
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Setup
xp=javax.xml.xpath.XPathFactory.newInstance.newXPath();

num_cal_vecs=str2double(xp.evaluate(...
    'count(calibration/calibrationVectorList/calibrationVector)',domnode));
betas = [];
for i = 1:num_cal_vecs
    betas = [betas str2num(xp.evaluate(...
        ['calibration/calibrationVectorList/calibrationVector[' num2str(i) ']/betaNought'],...
        domnode))];
end
last_az_time = xp.evaluate(...
        ['calibration/calibrationVectorList/calibrationVector[' num2str(num_cal_vecs) ']/azimuthTime'],...
        domnode);
if datenum(char(last_az_time),'yyyy-mm-ddTHH:MM:SS')<datenum(2015, 11, 25)
    warning('META2SICD_S1CAL:EARLY_DATE','Sentinel-1 radiometric data prior to November 25, 2015 might not be accurate.');    
end
if all(betas==betas(1))
    % Generally beta (and thus RCS) is constant across most S1 datasets
    % that we have seen, so we only handle the constant case here.
    output_meta.Radiometric.BetaZeroSFPoly = 1/betas(1)^2;
    % Other radiometric fields can be populated based off of this, but
    % this requires additional metadata, so that will be done outside this
    % function.
else
    output_meta = struct(); % No metadata to return if assumption not met
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////