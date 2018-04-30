function [ ccd_out, phase_out ] = ccdmem_angle( reference_image, match_image, corr_window_size, ang , angoutflag, NoiseVal, CCDMetric)
%CCDMEM_ANGLE - Copmputes CCD of two registered images, can input angle

if mod(ang,90) == 0
    flag = 0;
    if ang == 90 || ang == 270
        corr_window_size = fliplr(corr_window_size);
    end        
else
    flag = 1;
    [ny1,nx1] = size(reference_image);
    reference_image = imrotate(reference_image,ang);
    match_image = imrotate(match_image,ang);    
end

if ~exist('NoiseVal','var') || isempty(NoiseVal)
    NoiseVal = 0;
end

if ~exist('CCDMetric','var') || isempty(CCDMetric)
    CCDMetric = [];
end

if nargout>1
    if NoiseVal == 0
        [ccd_out,phase_out] = ccdmem(reference_image,match_image,corr_window_size);
    else
        if strcmpi(CCDMetric,'SCCM (Mitchell)')
            [ccd_out,phase_out] = SCCM(reference_image,match_image,corr_window_size(1),NoiseVal,NoiseVal);
        else
            [ccd_out,phase_out] = ccdnoisemem(reference_image,match_image,corr_window_size,NoiseVal,NoiseVal);
        end
    end
else
    if NoiseVal == 0
        ccd_out = ccdmem(reference_image,match_image,corr_window_size);
    else
        if strcmpi(CCDMetric,'SCCM (Mitchell)')
            ccd_out = SCCM(reference_image,match_image,corr_window_size(1),NoiseVal,NoiseVal);
        else
            ccd_out = ccdnoisemem(reference_image,match_image,corr_window_size,NoiseVal,NoiseVal);
        end
    end
end

if angoutflag==0   
    if flag
        ccd_out2 = imrotate(ccd_out,-1*ang);
        if nargout>1
            phase_out2 = imrotate(phase_out,-1*ang);
        end
        [ny2,nx2] = size(ccd_out2);
        clear ccd_out phase_out;
        xoffset = floor((nx2-nx1)/2);
        yoffset = floor((ny2-ny1)/2);
        ccd_out = ccd_out2(yoffset:(yoffset+ny1-1),xoffset:(xoffset+nx1-1));
        if nargout>1
            phase_out = phase_out2(yoffset:(yoffset+ny1-1),xoffset:(xoffset+nx1-1));
        end
    end       
else
    if ang == 90
        ccd_out = rot90(ccd_out,1);  
        if nargout>1
            phase_out = rot90(phase_out,1);
        end
    elseif ang == 270
        ccd_out = rot90(ccd_out,3);
        if nargout>1
            phase_out = rot90(phase_out,3);
        end
    end
end

end
