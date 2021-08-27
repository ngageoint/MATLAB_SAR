function [ output_meta ] = derived_sicd_fields( input_meta )
%DERIVED_SICD_FIELDS Computes derived fields from fundamental ones
%
% This function attempts to populate missing fields from a SICD metadata
% structure.  Using this function should allow one to more simply (with
% less replicated code) create SICDs from a number of different sources by
% defining only the fundamental fields and then calling this function to
% populate all of the derived fields.
% 
% There are two types of fields which are populated in this function:
% 1) DERIVED values: These fields can be computed exactly from other
% fields. SICD includes many redundant parameters for ease of access.  This
% function tries to see which core, fundamental fields are available and
% calculate as many derived fields from those as possible.  Example:
% SCPCOA.SCPTime must equal Grid.TimeCOAPoly(1).
% 2) DEFAULT values: These are fields which may not be given exactly, but
% for which we can make a reasonable guess or approximation based on the
% most common types of SAR processing.  In fact, some of these fields are
% so common that they are just assumed and not even explicitly given in
% other file formats. Population of these fields can be turned off through
% the SET_DEFAULT_VALUES variable since they are not absolutely known.
% Example: The PFA image plane normal is often the instantaneous slant
% plane at center of aperture.
%
% Within the code and comments, we attempt to label whether the value being
% populated is a DERIVED or DEFAULT value.  Note that if a field is already
% populated in the input metadata structure, this function will not
% overwrite it for either the DERIVED or DEFAULT cases.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

SET_DEFAULT_VALUES = true; % A default for using defaults...

output_meta=input_meta;

% Fields DERIVED from Grid parameters
for i = {'Row','Col'}
    row_column = i{1};
    if isfield(output_meta,'Grid')&&isfield(output_meta.Grid,row_column)
        % WgtFunct DERIVED from WgtType
        if isfield(output_meta.Grid.(row_column),'WgtType') && ...
                isfield(output_meta.Grid.(row_column).WgtType,'WindowName')&&...
                ~isfield(output_meta.Grid.(row_column),'WgtFunct') && ...
                ~any(strcmpi(output_meta.Grid.(row_column).WgtType.WindowName,{'UNIFORM','UNKNOWN'}))
            try % Will error if WgtFunct cannot be created from WgtType
                DEFAULT_WGT_SIZE = 512;
                WgtFunct = sicdweight2fun(output_meta.Grid.(row_column));
                output_meta.Grid.(row_column).WgtFunct = WgtFunct(DEFAULT_WGT_SIZE);
            end
        end
        if isfield(output_meta.Grid.(row_column),'WgtType') && ...
                isfield(output_meta.Grid.(row_column).WgtType,'WindowName')&&...
                strcmpi(output_meta.Grid.(row_column).WgtType.WindowName,'UNIFORM')
            broadening_factor = 2 * fzero(@(x) (sin(pi*x)/(pi*x)) - (1/sqrt(2)), .1); % 0.886
        elseif isfield(output_meta.Grid.(row_column),'WgtFunct')
            OVERSAMPLE = 1024;
            WgtFunct = output_meta.Grid.(row_column).WgtFunct;
            imp_resp = abs(fft(WgtFunct, round(numel(WgtFunct)*OVERSAMPLE))); % Oversampled response function
            imp_resp = imp_resp/sum(WgtFunct); % Normalize to unit peak
            ind = find(imp_resp<1/sqrt(2),1,'first')+[-1 -0]; % Samples surrounding half-power point
            ind = interp1(imp_resp(ind), ind, 1/sqrt(2)); % Linear interpolation to solve for half-power point
            broadening_factor = 2*(ind - 1)/OVERSAMPLE;
        end
        % Resolution can be DERIVED from bandwidth and weighting type
        if isfield(output_meta.Grid.(row_column),'ImpRespBW') && ...
                ~isfield(output_meta.Grid.(row_column),'ImpRespWid') && ...
                exist('broadening_factor','var')
            output_meta.Grid.(row_column).ImpRespWid=...
                broadening_factor/output_meta.Grid.(row_column).ImpRespBW;
        elseif isfield(output_meta.Grid.(row_column),'ImpRespWid') && ...
                ~isfield(output_meta.Grid.(row_column),'ImpRespBW') && ...
                exist('broadening_factor','var')
            output_meta.Grid.(row_column).ImpRespBW=...
                broadening_factor/output_meta.Grid.(row_column).ImpRespWid;
        end
        % DeltaK1/2 are APPROXIMATED from DeltaKCOAPoly
        if all(isfield(output_meta.Grid,{'Row','Col'}))&&...
                isfield(output_meta.Grid.(row_column),'DeltaKCOAPoly')&&...
                isfield(output_meta.Grid.(row_column),'ImpRespBW')&&...
                isfield(output_meta.Grid.(row_column),'SS')&&...
                ~any(isfield(output_meta.Grid.(row_column),{'DeltaK1','DeltaK2'}))
            % Here, we assume the min and max of DeltaKCOAPoly must be on
            % the vertices of the image, since it is smooth and monotonic
            % in most cases-- although in actuality this is not always the
            % case.  To be totally generic, we would have to search for an
            % interior min and max as well.
            if isfield(output_meta.ImageData,'ValidData') % Test vertices
                vertices = double([output_meta.ImageData.ValidData.Vertex.Col;...
                    output_meta.ImageData.ValidData.Vertex.Row]);
            else % Use edges of full image
                vertices = double([0 output_meta.ImageData.NumCols-1 ...
                    output_meta.ImageData.NumCols-1 0; ...
                    0 0 output_meta.ImageData.NumRows-1 ...
                    output_meta.ImageData.NumRows-1]);
            end
            if isfield(output_meta.Grid.(row_column),'DeltaKCOAPoly')
                min_dk = Inf; max_dk= -Inf;
            else
                min_dk = 0; max_dk= 0;
            end
            for j = 1:size(vertices,2)
                if isfield(output_meta.Grid.(row_column),'DeltaKCOAPoly')
                    currentDeltaK = sicd_polyval2d(output_meta.Grid.(row_column).DeltaKCOAPoly,...
                        vertices(1,j),vertices(2,j),output_meta);
                    min_dk = min(min_dk, currentDeltaK);
                    max_dk = max(max_dk, currentDeltaK);
                end
            end
            min_dk = min_dk - (output_meta.Grid.(row_column).ImpRespBW/2);
            max_dk = max_dk + (output_meta.Grid.(row_column).ImpRespBW/2);
            % Wrapped spectrum
            if (min_dk < -(1/output_meta.Grid.(row_column).SS)/2) || ...
                    (max_dk > (1/output_meta.Grid.(row_column).SS)/2)
                min_dk = -(1/output_meta.Grid.(row_column).SS)/2;
                max_dk = -min_dk;
            end
            output_meta.Grid.(row_column).DeltaK1=min_dk;
            output_meta.Grid.(row_column).DeltaK2=max_dk;
        end
    end
end

% SCPTime can always be DERIVED from Grid.TimeCOAPoly
if (~isfield(output_meta,'SCPCOA')||~isfield(output_meta.SCPCOA,'SCPTime')) && ...
        isfield(output_meta,'Grid')&&isfield(output_meta.Grid,'TimeCOAPoly')
    output_meta.SCPCOA.SCPTime = output_meta.Grid.TimeCOAPoly(1,1);
% and sometimes Grid.TimeCOAPoly can be DERIVED from SCPTime
elseif (~isfield(output_meta,'Grid')||~isfield(output_meta.Grid,'TimeCOAPoly')) && ...
        isfield(output_meta,'SCPCOA')&&isfield(output_meta.SCPCOA,'SCPTime')&&...
        isfield(output_meta,'CollectionInfo')&&...
        isfield(output_meta.CollectionInfo,'RadarMode')&&...
        isfield(output_meta.CollectionInfo.RadarMode,'ModeType')&&...
        strcmpi(output_meta.CollectionInfo.RadarMode.ModeType,'SPOTLIGHT')
    output_meta.Grid.TimeCOAPoly=output_meta.SCPCOA.SCPTime;
end

% ARP Pos/Vel/ACC fields can be DERIVED from ARPPoly and SCPTime
if isfield(output_meta,'Position')&&isfield(output_meta.Position,'ARPPoly')&&...
        isfield(output_meta,'SCPCOA')&&isfield(output_meta.SCPCOA,'SCPTime')
    % ARPPoly should always be a vertical (column) vector
    pos_coefs=[output_meta.Position.ARPPoly.X(end:-1:1)...
        output_meta.Position.ARPPoly.Y(end:-1:1)...
        output_meta.Position.ARPPoly.Z(end:-1:1)];
    arp_pos_ecf=[polyval(pos_coefs(:,1),output_meta.SCPCOA.SCPTime)...
        polyval(pos_coefs(:,2),output_meta.SCPCOA.SCPTime)...
        polyval(pos_coefs(:,3),output_meta.SCPCOA.SCPTime)];
    if ~(isfield(output_meta,'SCPCOA')&&isfield(output_meta.SCPCOA,'ARPPos'))
        output_meta.SCPCOA.ARPPos.X=arp_pos_ecf(1);
        output_meta.SCPCOA.ARPPos.Y=arp_pos_ecf(2);
        output_meta.SCPCOA.ARPPos.Z=arp_pos_ecf(3);
    end
    % Velocity is derivative of position
    vel_coefs=pos_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-1):-1:1)',[1 3]);
    arp_vel_ecf=[polyval(vel_coefs(:,1), output_meta.SCPCOA.SCPTime)...
        polyval(vel_coefs(:,2), output_meta.SCPCOA.SCPTime)...
        polyval(vel_coefs(:,3), output_meta.SCPCOA.SCPTime)];
    if ~(isfield(output_meta,'SCPCOA')&&isfield(output_meta.SCPCOA,'ARPVel'))
        output_meta.SCPCOA.ARPVel.X=arp_vel_ecf(1);
        output_meta.SCPCOA.ARPVel.Y=arp_vel_ecf(2);
        output_meta.SCPCOA.ARPVel.Z=arp_vel_ecf(3);
    end
    if ~(isfield(output_meta,'SCPCOA')&&isfield(output_meta.SCPCOA,'ARPAcc'))
        % Acceleration is derivative of velocity
        acc_coefs=vel_coefs(1:end-1,:).*repmat(((size(pos_coefs,1)-2):-1:1)',[1 3]);
        output_meta.SCPCOA.ARPAcc.X=...
            polyval(acc_coefs(:,1),output_meta.SCPCOA.SCPTime);
        output_meta.SCPCOA.ARPAcc.Y=...
            polyval(acc_coefs(:,2),output_meta.SCPCOA.SCPTime);
        output_meta.SCPCOA.ARPAcc.Z=...
            polyval(acc_coefs(:,3),output_meta.SCPCOA.SCPTime);
    end
end
% A simple ARPPoly can be DERIVED from SCPCOA Pos/Vel/Acc if that was
% all that was defined.
if isfield(output_meta,'SCPCOA') && ...
        all(isfield(output_meta.SCPCOA,{'ARPPos','ARPVel','SCPTime'})) && ...
        (~isfield(output_meta,'Position')||~isfield(output_meta.Position,'ARPPoly'))
    if ~isfield(output_meta.SCPCOA,'ARPAcc')
        output_meta.SCPCOA.ARPAcc = struct('X',0,'Y',0,'Z',0);
    end
    for i = {'X','Y','Z'}
        output_meta.Position.ARPPoly.(i{1}) = [...
            ... % Constant
            output_meta.SCPCOA.ARPPos.(i{1}) - ...
            (output_meta.SCPCOA.ARPVel.(i{1}) * output_meta.SCPCOA.SCPTime) + ...
            ((output_meta.SCPCOA.ARPAcc.(i{1})/2) * (output_meta.SCPCOA.SCPTime^2)); ...
            ... % Linear
            output_meta.SCPCOA.ARPVel.(i{1}) - ...
            (output_meta.SCPCOA.ARPAcc.(i{1}) * output_meta.SCPCOA.SCPTime); ...
            ... % Quadratic
            output_meta.SCPCOA.ARPAcc.(i{1})/2];
    end
end

% Transmit bandwidth
if isfield(output_meta,'RadarCollection') && ...
        isfield(output_meta.RadarCollection,'Waveform') && ...
        isfield(output_meta.RadarCollection.Waveform,'WFParameters')
    % DERIVED: Redundant WFParameters fields
    temp_wf = output_meta.RadarCollection.Waveform; % Adding fields to output_meta will add empty values in all WFParameters after the first
    for i = 1:numel(output_meta.RadarCollection.Waveform.WFParameters)
        wfp = temp_wf.WFParameters(i); % To shorten notation
        if isfield(wfp, 'RcvDemodType') && strcmp(wfp.RcvDemodType, 'CHIRP') && ~isfield(wfp, 'RcvFMRate') 
            output_meta.RadarCollection.Waveform.WFParameters(i).RcvFMRate = 0;
        end
        if isfield(wfp, 'RcvFMRate') && (wfp.RcvFMRate==0) && ~isfield(wfp, 'RcvDemodType') 
            output_meta.RadarCollection.Waveform.WFParameters(i).RcvDemodType = 'CHIRP';
        end
        if ~isfield(wfp, 'TxRFBandwidth') && isfield(wfp, 'TxPulseLength') && isfield(wfp, 'TxFMRate')
            output_meta.RadarCollection.Waveform.WFParameters(i).TxRFBandwidth = ...
                wfp.TxPulseLength * wfp.TxFMRate;
        end
        if isfield(wfp, 'TxRFBandwidth') && ~isfield(wfp, 'TxPulseLength') && isfield(wfp, 'TxFMRate')
            output_meta.RadarCollection.Waveform.WFParameters(i).TxPulseLength = ...
                wfp.TxRFBandwidth / wfp.TxFMRate;
        end
        if isfield(wfp, 'TxRFBandwidth') && isfield(wfp, 'TxPulseLength') && ~isfield(wfp, 'TxFMRate')
            output_meta.RadarCollection.Waveform.WFParameters(i).TxFMRate = ...
                wfp.TxRFBandwidth / wfp.TxPulseLength;
        end
    end
    % DERIVED: These values should be equal.
    if ~isfield(output_meta.RadarCollection,'TxFrequency') || ...
            ~isfield(output_meta.RadarCollection.TxFrequency,'Min')
        output_meta.RadarCollection.TxFrequency.Min = ...
            min([output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart]);
    end
    if ~isfield(output_meta.RadarCollection,'TxFrequency') || ...
            ~isfield(output_meta.RadarCollection.TxFrequency,'Max')
        output_meta.RadarCollection.TxFrequency.Max = ...
            max([output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart] + ...
            [output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth]);
    end
end
if isfield(output_meta,'RadarCollection') && ...
        isfield(output_meta.RadarCollection,'TxFrequency') && ...
        all(isfield(output_meta.RadarCollection.TxFrequency,{'Min','Max'}))
    % DEFAULT: We often assume that all transmitted bandwidth was
    % processed, if given no other information.
    if SET_DEFAULT_VALUES
        if ~isfield(output_meta,'ImageFormation') || ...
                ~isfield(output_meta.ImageFormation,'TxFrequencyProc') || ...
                ~isfield(output_meta.ImageFormation.TxFrequencyProc,'MinProc')
            output_meta.ImageFormation.TxFrequencyProc.MinProc = output_meta.RadarCollection.TxFrequency.Min;
        end
        if ~isfield(output_meta,'ImageFormation') || ...
                ~isfield(output_meta.ImageFormation,'TxFrequencyProc') || ...
                ~isfield(output_meta.ImageFormation.TxFrequencyProc,'MaxProc')
            output_meta.ImageFormation.TxFrequencyProc.MaxProc = output_meta.RadarCollection.TxFrequency.Max;
        end
    end
    
    % DERIVED: These values should be equal.
    if isfield(output_meta.RadarCollection,'Waveform') && ...
            isfield(output_meta.RadarCollection.Waveform,'WFParameters') && ...
            numel(output_meta.RadarCollection.Waveform.WFParameters)==1
        if ~isfield(output_meta.RadarCollection.Waveform.WFParameters,'TxFreqStart')
            output_meta.RadarCollection.Waveform.WFParameters.TxFreqStart =...
                output_meta.RadarCollection.TxFrequency.Min;
        end
        if ~isfield(output_meta.RadarCollection.Waveform.WFParameters,'TxRFBandwidth')
            output_meta.RadarCollection.Waveform.WFParameters.TxRFBandwidth =...
                output_meta.RadarCollection.TxFrequency.Max -...
                output_meta.RadarCollection.TxFrequency.Min;
        end
    end
end

% We might use center processed frequency later
if isfield(output_meta,'ImageFormation') && ...
        isfield(output_meta.ImageFormation,'TxFrequencyProc') && ...
        all(isfield(output_meta.ImageFormation.TxFrequencyProc,{'MinProc','MaxProc'})) && ...
        (~isfield(output_meta,'RadarCollection') || ...
        ~isfield(output_meta.RadarCollection,'RefFreqIndex') || ...
        (output_meta.RadarCollection.RefFreqIndex==0))
    fc = (output_meta.ImageFormation.TxFrequencyProc.MinProc + ...
        output_meta.ImageFormation.TxFrequencyProc.MaxProc)/2;
end

% DERIVED: GeoData.SCP
if isfield(output_meta, 'GeoData') && isfield(output_meta.GeoData,'SCP') && ...
        isfield(output_meta.GeoData.SCP,'ECF') && ~isfield(output_meta.GeoData.SCP,'LLH')
    llh=ecf_to_geodetic([output_meta.GeoData.SCP.ECF.X ...
        output_meta.GeoData.SCP.ECF.Y output_meta.GeoData.SCP.ECF.Z]);
    output_meta.GeoData.SCP.LLH.Lat=llh(1);
    output_meta.GeoData.SCP.LLH.Lon=llh(2);
    output_meta.GeoData.SCP.LLH.HAE=llh(3);
end
if isfield(output_meta, 'GeoData') && isfield(output_meta.GeoData,'SCP') && ...
        isfield(output_meta.GeoData.SCP,'LLH') && ~isfield(output_meta.GeoData.SCP,'ECF')
    ecf=geodetic_to_ecf([output_meta.GeoData.SCP.LLH.Lat ...
        output_meta.GeoData.SCP.LLH.Lon output_meta.GeoData.SCP.LLH.HAE]);
    output_meta.GeoData.SCP.ECF.X=ecf(1);
    output_meta.GeoData.SCP.ECF.Y=ecf(2);
    output_meta.GeoData.SCP.ECF.Z=ecf(3);
end

% Many fields (particularly in SCPCOA) can be DERIVED from ARPPos, ARPVel and SCP
if isfield(output_meta,'SCPCOA')&&all(isfield(output_meta.SCPCOA,{'ARPPos','ARPVel'}))&&...
        isfield(output_meta,'GeoData')&&isfield(output_meta.GeoData,'SCP')&&...
        isfield(output_meta.GeoData.SCP,'ECF')
    % GroundRange computation is particularly sensitive to precision, so we
    % wrap all of these in double, just in case they were passed in as
    % single.  The values themselves don't need double precision, but the
    % computation does.
    SCP=double([output_meta.GeoData.SCP.ECF.X output_meta.GeoData.SCP.ECF.Y output_meta.GeoData.SCP.ECF.Z]);
    ARP=double([output_meta.SCPCOA.ARPPos.X output_meta.SCPCOA.ARPPos.Y output_meta.SCPCOA.ARPPos.Z]);
    ARP_v=double([output_meta.SCPCOA.ARPVel.X output_meta.SCPCOA.ARPVel.Y output_meta.SCPCOA.ARPVel.Z]);
    uLOS = (SCP - ARP).'/norm(SCP - ARP);
    left = cross(ARP/norm(ARP),ARP_v/norm(ARP_v));
    look = sign(left * uLOS);
    if ~isfield(output_meta.SCPCOA,'SideOfTrack')
        if look<0
            output_meta.SCPCOA.SideOfTrack='R';
        else
            output_meta.SCPCOA.SideOfTrack='L';
        end
    end
    if ~isfield(output_meta.SCPCOA,'SlantRange')
        output_meta.SCPCOA.SlantRange = norm(SCP - ARP);
    end
    if ~isfield(output_meta.SCPCOA,'GroundRange')
        output_meta.SCPCOA.GroundRange = ...
            norm(SCP) * acos(ARP * SCP.' / (norm(SCP) * norm(ARP)));
    end
    if ~isfield(output_meta.SCPCOA,'DopplerConeAng')
         % Doppler Cone Angle is angle of slant range vector from velocity vector
        output_meta.SCPCOA.DopplerConeAng = ...
            acosd((ARP_v / norm(ARP_v)) * uLOS);
    end
    % Earth Tangent Plane (ETP) at the SCP is the plane tangent to the
    % surface of constant height above the WGS 84 ellipsoid (HAE) that
    % contains the SCP. The ETP is an approximation to the ground plane at
    % the SCP.
    ETP=wgs_84_norm(SCP).';
    if ~isfield(output_meta.SCPCOA,'GrazeAng')
        % Angle between ground plane and line-of-site vector
        output_meta.SCPCOA.GrazeAng = asind(ETP * (-uLOS));
    end
    if ~isfield(output_meta.SCPCOA,'IncidenceAng')
        % Angle between ground plane normal and line-of-site vector
        % output_meta.SCPCOA.IncidenceAng = acosd(ETP * (-uLOS));
        output_meta.SCPCOA.IncidenceAng = 90 - output_meta.SCPCOA.GrazeAng;
    end
    spn=look*cross(ARP_v,uLOS); spn=spn/norm(spn); % Instantaneous slant plane unit normal at COA (also called uSPZ in SICD spec)
    uGPX = (-uLOS) - (ETP * (-uLOS)) * ETP.'; % Project range vector (from SCP toward ARP) onto ground plane
    uGPX = uGPX/norm(uGPX);
    if ~isfield(output_meta.SCPCOA,'TwistAng')
        % 1) Equations from SICD spec:
        uGPY = cross(ETP,uGPX);
        % Angle from the +GPY axis and to the +SPY axis in the plane of incidence
        output_meta.SCPCOA.TwistAng = -asind(uGPY*spn.');
        % 2) Another implementation (seems to turn out exactly the same):
        % output_meta.SCPCOA.TwistAng = asind(cross(ETP, spn) * (-uLOS) / norm(cross((-uLOS), ETP)));
    end
    if ~isfield(output_meta.SCPCOA,'SlopeAng')
        output_meta.SCPCOA.SlopeAng = acosd(ETP * spn.'); % Angle between slant and ground planes
    end
    north_ground = [ 0; 0; 1 ] - (ETP * [ 0; 0; 1 ]) * ETP.'; % Project north onto ground plane
    uNORTH = north_ground/norm(north_ground); % Unit vector in ground plane in north direction
    uEAST = cross(uNORTH, ETP); % Unit vector in ground plane in east direction
    if ~isfield(output_meta.SCPCOA,'AzimAng')
        az_north = uGPX.' * uNORTH; % Component of ground-projected range vector in north direction
        az_east = uGPX.' * uEAST.'; % Component of ground-projected range vector in east direction
        output_meta.SCPCOA.AzimAng = atan2( az_east, az_north);
        output_meta.SCPCOA.AzimAng = mod(output_meta.SCPCOA.AzimAng*180/pi, 360); % Assure in [0,360], not [-pi,pi]
    end
    if ~isfield(output_meta.SCPCOA,'LayoverAng')
        layover_ground = (ETP - (1 / (ETP * spn.')) * spn).'; % Layover direction in ground plane
        lo_north = layover_ground.' * uNORTH; % Component of layover in north direction
        lo_east = layover_ground.' * uEAST.'; % Component of layover in east direction
        output_meta.SCPCOA.LayoverAng = atan2( lo_east, lo_north);
        output_meta.SCPCOA.LayoverAng = mod(output_meta.SCPCOA.LayoverAng*180/pi, 360); % Assure in [0,360], not [-pi,pi]
    end

    % Compute IFP specific parameters (including Row/Col.UVectECF) here
    if isfield(output_meta,'ImageFormation')&&...
            isfield(output_meta.ImageFormation,'ImageFormAlgo')&&...
            ischar(output_meta.ImageFormation.ImageFormAlgo)
        switch upper(output_meta.ImageFormation.ImageFormAlgo)
            case 'RGAZCOMP'
                % In theory, we could even derive Grid.TimeCOAPoly for the
                % RGAZCOMP case if IPPPoly was include, since it must just
                % be the time computed for the vector index: v_coa = (1/2)
                % * (v_ps + v_pe)
                % DERIVED: RGAZCOMP image formation must result in a SLANT, RGAZIM grid
                if ~isfield(output_meta,'Grid') || ~isfield(output_meta.Grid,'ImagePlane')
                    output_meta.Grid.ImagePlane = 'SLANT';
                end
                if ~isfield(output_meta,'Grid') || ~isfield(output_meta.Grid,'Type')
                    output_meta.Grid.Type = 'RGAZIM';
                end
                % DERIVED: RgAzComp.AzSF
                if ~isfield(output_meta,'RgAzComp') || ...
                        ~isfield(output_meta.RgAzComp,'AzSF')
                    output_meta.RgAzComp.AzSF = ...
                        -look * sind(output_meta.SCPCOA.DopplerConeAng)/...
                        output_meta.SCPCOA.SlantRange;
                end
                % DERIVED: RgAzComp.KazPoly
                if isfield(output_meta,'Timeline') && ...
                        isfield(output_meta.Timeline,'IPP') && ...
                        isfield(output_meta.Timeline.IPP,'Set') && ...
                        numel(output_meta.Timeline.IPP.Set)==1 && ...
                        isfield(output_meta.Timeline.IPP.Set,'IPPPoly') && ...
                        isfield(output_meta,'Grid') && ...
                        isfield(output_meta.Grid,'Row') && ...
                        isfield(output_meta.Grid.Row,'KCtr') && ...
                        ~isfield(output_meta.RgAzComp,'KazPoly')
                    krg_coa = output_meta.Grid.Row.KCtr;
                    if isfield(output_meta.Grid.Row,'DeltaKCOAPoly')
                        krg_coa = krg_coa + output_meta.Grid.Row.DeltaKCOAPoly;
                    end
                    st_rate_coa = polyval(polyder(output_meta.Timeline.IPP.Set.IPPPoly(end:-1:1)),output_meta.SCPCOA.SCPTime);
                    delta_kaz_per_delta_v = look * krg_coa * ... % Scale factor described in SICD spec
                        (norm(ARP_v) * sind(output_meta.SCPCOA.DopplerConeAng) / output_meta.SCPCOA.SlantRange) / ...
                        st_rate_coa;
                    output_meta.RgAzComp.KazPoly = delta_kaz_per_delta_v * ...
                        output_meta.Timeline.IPP.Set.IPPPoly;
                end
                % DERIVED: UVectECF
                if (~isfield(output_meta,'Grid') || ...
                        ((~isfield(output_meta.Grid,'Row') || ~isfield(output_meta.Grid.Row,'UVectECF')) && ...
                        (~isfield(output_meta.Grid,'Col') || ~isfield(output_meta.Grid.Col,'UVectECF'))))
                    output_meta.Grid.Row.UVectECF = ...
                        cell2struct(num2cell(uLOS),{'X','Y','Z'});
                    output_meta.Grid.Col.UVectECF = ...
                        cell2struct(num2cell(cross(spn(:),uLOS)),{'X','Y','Z'});
                end
                % DERIVED: KCtr/DeltaKCOAPoly
                % In SICD, if the optional DeltaKCOAPoly field is omitted,
                % it is assumed to be zero. If the creator of the partial
                % SICD metadata just forgot it, or didn't know it, rather
                % than leaving the field off as an explicit declaration of
                % a zero value, the KCtr computation will be wrong if the
                % DFT was not "centered" (s_0 = s_coa and v_0 = v_coa in
                % the terminology of the SICD spec).
                if exist('fc','var')
                    if ~isfield(output_meta,'Grid') || ...
                            ~isfield(output_meta.Grid,'Row') || ...
                            ~isfield(output_meta.Grid.Row,'KCtr')
                        if isfield(output_meta,'Grid') && ...
                            isfield(output_meta.Grid,'Row') && ...
                            isfield(output_meta.Grid.Row,'DeltaKCOAPoly')
                            % DeltaKCOAPoly populated, but not KCtr (would be odd)
                            output_meta.Grid.Row.KCtr = (fc * (2/SPEED_OF_LIGHT)) - ...
                                output_meta.Grid.Row.DeltaKCOAPoly(1);
                        else % Neither KCtr or DeltaKCOAPoly populated
                            output_meta.Grid.Row.KCtr = fc * (2/SPEED_OF_LIGHT);
                            % DeltaKCOAPoly assumed to be zero
                        end
                    elseif ~isfield(output_meta.Grid.Row,'DeltaKCOAPoly')
                        % KCtr populated, but not DeltaKCOAPoly
                        output_meta.Grid.Row.DeltaKCOAPoly = (fc * (2/SPEED_OF_LIGHT)) - ...
                            output_meta.Grid.Row.KCtr;
                    end
                end
                if ~isfield(output_meta,'Grid') || ...
                        ~isfield(output_meta.Grid,'Col') || ...
                        ~isfield(output_meta.Grid.Col,'KCtr')
                    output_meta.Grid.Col.KCtr = 0;
                    if isfield(output_meta.Grid.Col,'DeltaKCOAPoly')
                        % DeltaKCOAPoly populated, but not KCtr (would be odd)
                        output_meta.Grid.Col.KCtr = -output_meta.Grid.Col.DeltaKCOAPoly(1);
                    else % Neither KCtr or DeltaKCOAPoly populated
                        % DeltaKCOAPoly assumed to be zero
                    end
                elseif ~isfield(output_meta.Grid.Col,'DeltaKCOAPoly')
                    % KCtr populated, but not DeltaKCOAPoly
                    output_meta.Grid.Col.DeltaKCOAPoly = -output_meta.Grid.Col.KCtr;
                end
            case 'PFA'
                % DEFAULT: RGAZIM grid is the natural result of PFA
                if SET_DEFAULT_VALUES && (~isfield(output_meta,'Grid') || ~isfield(output_meta.Grid,'Type'))
                    output_meta.Grid.Type = 'RGAZIM';
                end
                % Reasonable DEFAULT guesses for PFA parameters IPN, FPN,
                % and PolarAngRefTime
                if SET_DEFAULT_VALUES && (~isfield(output_meta,'PFA') || ~isfield(output_meta.PFA,'IPN'))
                    if isfield(output_meta,'Grid')&&isfield(output_meta.Grid,'ImagePlane')
                        if strcmpi(output_meta.Grid.ImagePlane,'SLANT')
                            % Instantaneous slant plane at center of aperture
                            output_meta.PFA.IPN=cell2struct(num2cell(spn.'),{'X','Y','Z'});
                        elseif strcmpi(output_meta.Grid.ImagePlane,'GROUND')
                            output_meta.PFA.IPN=cell2struct(num2cell(ETP.'),{'X','Y','Z'});
                        end
                    else % Guess slant plane (the most common IPN) if not specified
                        output_meta.PFA.IPN=cell2struct(num2cell(spn.'),{'X','Y','Z'});
                    end
                end
                if SET_DEFAULT_VALUES && (~isfield(output_meta,'PFA') || ~isfield(output_meta.PFA,'FPN'))
                    output_meta.PFA.FPN=cell2struct(num2cell(ETP.'),{'X','Y','Z'});
                end
                if isfield(output_meta,'Position') && ...
                        isfield(output_meta.Position,'ARPPoly') && ...
                        isfield(output_meta,'PFA') && ...
                        isfield(output_meta.PFA,'PolarAngRefTime') % Compute exactly if possible
                    pol_ref_pos = [polyval(output_meta.Position.ARPPoly.X(end:-1:1), output_meta.PFA.PolarAngRefTime) ...
                        polyval(output_meta.Position.ARPPoly.Y(end:-1:1), output_meta.PFA.PolarAngRefTime) ...
                        polyval(output_meta.Position.ARPPoly.Z(end:-1:1), output_meta.PFA.PolarAngRefTime)];
                elseif SET_DEFAULT_VALUES % DEFAULT: Otherwise guess PolarAngRefTime = SCPTime
                    pol_ref_pos = ARP;
                    if isfield(output_meta,'SCPCOA') && isfield(output_meta.SCPCOA,'SCPTime')
                        output_meta.PFA.PolarAngRefTime = output_meta.SCPCOA.SCPTime;
                    end
                end
                % A following fields can be DERIVED from those above
                if isfield(output_meta,'Timeline') && isfield(output_meta, 'PFA') && ...
                        isfield(output_meta.Timeline,'CollectDuration') && ...
                        ~any(isfield(output_meta.PFA,{'PolarAngPoly','SpatialFreqSFPoly'})) && ...
                        isfield(output_meta.Position,'ARPPoly')
                    % Compute values we will need for pfa_polar_coords.
                    % Sample the collection period evenly.  Don't need too
                    % many sample points since the ARPPoly is typically a
                    % fairly low order polynomial (<=5).
                    times = linspace(0, output_meta.Timeline.CollectDuration, 15).';
                    pos = [polyval(output_meta.Position.ARPPoly.X(end:-1:1), times) ...
                        polyval(output_meta.Position.ARPPoly.Y(end:-1:1), times) ...
                        polyval(output_meta.Position.ARPPoly.Z(end:-1:1), times)];
                    [k_a, k_sf] = pfa_polar_coords(pos, SCP, pol_ref_pos, ...
                        [output_meta.PFA.IPN.X output_meta.PFA.IPN.Y output_meta.PFA.IPN.Z], ...
                        [output_meta.PFA.FPN.X output_meta.PFA.FPN.Y output_meta.PFA.FPN.Z]);
                    old_state = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
                    output_meta.PFA.PolarAngPoly = fliplr( polyfit(times, k_a, 5) ).';
                    output_meta.PFA.SpatialFreqSFPoly = fliplr( polyfit(k_a, k_sf, 5) ).';
                    warning(old_state);
                end

                if isfield(output_meta,'PFA') && ...
                        isfield(output_meta.PFA,'IPN') && ...
                        isfield(output_meta.PFA,'FPN') && ...
                        (~isfield(output_meta,'Grid') || ...
                        ((~isfield(output_meta.Grid,'Row') || ~isfield(output_meta.Grid.Row,'UVectECF')) && ...
                        (~isfield(output_meta.Grid,'Col') || ~isfield(output_meta.Grid.Col,'UVectECF'))))
                    ipn = [output_meta.PFA.IPN.X output_meta.PFA.IPN.Y output_meta.PFA.IPN.Z];
                    fpn = [output_meta.PFA.FPN.X output_meta.PFA.FPN.Y output_meta.PFA.FPN.Z];
                    % Row.UVect should be the range vector at zero polar
                    % angle time projected into IPN
                    % Projection of a point along a given direction to a
                    % plane is just the intersection of the line defined by
                    % that point (l0) and direction (l) and the plane
                    % defined by a point in the plane (p0) and the normal
                    % (p):
                    % l0 + ((l0 - p0).p/(l.p))*l
                    % where . represents the dot product.
                    d = (bsxfun(@minus, SCP, pol_ref_pos) * ipn(:)) ./...
                        (fpn * ipn(:)); % Distance from point to plane in line_direction
                    ref_pos_ipn = pol_ref_pos + (d * fpn);
                    uRG = (SCP - ref_pos_ipn).';
                    uRG = uRG/norm(uRG);
                    output_meta.Grid.Row.UVectECF = ...
                        cell2struct(num2cell(uRG),{'X','Y','Z'});
                    output_meta.Grid.Col.UVectECF = ...
                        cell2struct(num2cell(cross(ipn,uRG).'),{'X','Y','Z'});
                end
                % DEFAULT value. Almost always zero for PFA
                if SET_DEFAULT_VALUES && (~isfield(output_meta,'Grid') || ...
                        ~isfield(output_meta.Grid,'Col') || ...
                        ~isfield(output_meta.Grid.Col,'KCtr'))
                    output_meta.Grid.Col.KCtr = 0;
                    % Sometimes set to a nonzero (PFA.Kaz1 + PFA.Kaz2)/2
                end
                if SET_DEFAULT_VALUES && (~isfield(output_meta,'Grid') || ...
                        ~isfield(output_meta.Grid,'Row') || ...
                        ~isfield(output_meta.Grid.Row,'KCtr'))
                    if isfield(output_meta,'PFA') && ...
                            all(isfield(output_meta.PFA,{'Krg1','Krg2'}))
                        % DEFAULT: The most reasonable way to compute this
                        output_meta.Grid.Row.KCtr = (output_meta.PFA.Krg1 + output_meta.PFA.Krg2)/2;
                    elseif exist('fc','var')
                        % APPROXIMATION: This may not be quite right, due
                        % to rectangular inscription loss in PFA, but it
                        % should be close.
                        output_meta.Grid.Row.KCtr = fc * (2/SPEED_OF_LIGHT) * ...
                            output_meta.PFA.SpatialFreqSFPoly(1);
                    end
                end
            case 'RMA'
                if isfield(output_meta,'RMA') && isfield(output_meta.RMA,'ImageType')
                    rmatype = output_meta.RMA.ImageType;
                    % RMAT/RMCR cases
                    if any(strcmpi(rmatype,{'RMAT','RMCR'}))
                        if SET_DEFAULT_VALUES
                            if ~isfield(output_meta,'Grid') || ...
                                    ~isfield(output_meta.Grid,'ImagePlane')
                                output_meta.Grid.ImagePlane = 'SLANT';
                            end
                            if ~isfield(output_meta,'Grid') || ~isfield(output_meta.Grid,'Type')
                                % DEFAULT: XCTYAT grid is the natural result of RMA/RMAT
                                if strcmpi(rmatype,'RMAT')
                                    output_meta.Grid.Type = 'XCTYAT';
                                % DEFAULT: XRGYCR grid is the natural result of RMA/RMCR
                                elseif strcmpi(rmatype,'RMCR')
                                    output_meta.Grid.Type = 'XRGYCR';
                                end
                            end
                            % DEFAULT: Set PosRef/VelRef to SCPCOA Pos/Vel
                            if ~isfield(output_meta.RMA,rmatype) || ...
                                    ~isfield(output_meta.RMA.(rmatype),'PosRef')
                                output_meta.RMA.(rmatype).PosRef = output_meta.SCPCOA.ARPPos;
                            end
                            if ~isfield(output_meta.RMA,rmatype) || ...
                                    ~isfield(output_meta.RMA.(rmatype),'VelRef')
                                output_meta.RMA.(rmatype).VelRef = output_meta.SCPCOA.ARPVel;
                            end
                            % DEFAULT: RMAT/RMCR Row/Col.KCtr
                            if exist('fc','var')
                                k_f_c = fc * (2/SPEED_OF_LIGHT);
                                if strcmpi(rmatype,'RMAT')
                                    if ~isfield(output_meta,'Grid') || ...
                                            ~isfield(output_meta.Grid,'Row') || ...
                                            ~isfield(output_meta.Grid.Row,'KCtr')
                                        output_meta.Grid.Row.KCtr = k_f_c * ...
                                            sind(output_meta.RMA.RMAT.DopConeAngRef);
                                    end
                                    if ~isfield(output_meta,'Grid') || ...
                                            ~isfield(output_meta.Grid,'Col') || ...
                                            ~isfield(output_meta.Grid.Col,'KCtr')
                                        output_meta.Grid.Col.KCtr = k_f_c * ...
                                            cosd(output_meta.RMA.RMAT.DopConeAngRef);
                                    end
                                elseif strcmpi(rmatype,'RMCR')
                                    if ~isfield(output_meta,'Grid') || ...
                                            ~isfield(output_meta.Grid,'Row') || ...
                                            ~isfield(output_meta.Grid.Row,'KCtr')
                                        output_meta.Grid.Row.KCtr=k_f_c;
                                    end
                                    if ~isfield(output_meta,'Grid') || ...
                                            ~isfield(output_meta.Grid,'Col') || ...
                                            ~isfield(output_meta.Grid.Col,'KCtr')
                                        output_meta.Grid.Col.KCtr=0;
                                    end
                                end
                            end
                        end
                        if isfield(output_meta.RMA,rmatype) && ...
                                all(isfield(output_meta.RMA.(rmatype),{'PosRef','VelRef'}))
                            PosRef=[output_meta.RMA.(rmatype).PosRef.X ...
                                output_meta.RMA.(rmatype).PosRef.Y output_meta.RMA.(rmatype).PosRef.Z];
                            VelRef = [output_meta.RMA.(rmatype).VelRef.X ...
                                output_meta.RMA.(rmatype).VelRef.Y output_meta.RMA.(rmatype).VelRef.Z];
                            uLOS = (SCP - PosRef)/norm(SCP - PosRef); % Range unit vector
                            left = cross(PosRef/norm(PosRef),VelRef/norm(VelRef));
                            look = sign(left * uLOS');
                            % DCA is a DERIVED field
                            if ~isfield(output_meta.RMA.(rmatype),'DopConeAngRef')
                                output_meta.RMA.(rmatype).DopConeAngRef=acosd((VelRef/norm(VelRef))*uLOS.'); %
                            end
                            % Row/Col.UVectECF are DERIVED fields
                            if (~isfield(output_meta,'Grid') || ...
                                    ((~isfield(output_meta.Grid,'Row') || ...
                                    ~isfield(output_meta.Grid.Row,'UVectECF')) && ...
                                    (~isfield(output_meta.Grid,'Col') || ...
                                    ~isfield(output_meta.Grid.Col,'UVectECF'))))
                                if strcmp(rmatype,'RMAT')
                                    uYAT = -look*VelRef/norm(VelRef); % Along track unit vector
                                    spn = cross(uLOS,uYAT); spn=spn/norm(spn); % Reference slant plane normal
                                    uXCT = cross(uYAT, spn); % Cross track unit vector
                                    output_meta.Grid.Row.UVectECF = ...
                                        cell2struct(num2cell(uXCT),{'X','Y','Z'});
                                    output_meta.Grid.Col.UVectECF = ...
                                        cell2struct(num2cell(uYAT),{'X','Y','Z'});
                                elseif strcmp(rmatype,'RMCR')
                                    uXRG = uLOS; % Range unit vector
                                    spn = look * cross(VelRef / norm(VelRef), uXRG);
                                    spn=spn/norm(spn); % Reference slant plane normal
                                    uYCR = cross(spn, uXRG); % Cross range unit vector
                                    output_meta.Grid.Row.UVectECF = ...
                                        cell2struct(num2cell(uXRG(:)),{'X','Y','Z'});
                                    output_meta.Grid.Col.UVectECF = ...
                                        cell2struct(num2cell(uYCR(:)),{'X','Y','Z'});
                                end
                            end
                        end
                    % INCA
                    elseif strcmpi(output_meta.RMA.ImageType,'INCA') && ...
                            isfield(output_meta.RMA,'INCA')
                        % DEFAULT: RGZERO grid is the natural result of RMA/INCA
                        if ~isfield(output_meta,'Grid') || ~isfield(output_meta.Grid,'Type')
                            output_meta.Grid.Type = 'RGZERO';
                        end
                        if isfield(output_meta.RMA.INCA,'TimeCAPoly') && ...
                                isfield(output_meta,'Position') && ...
                                isfield(output_meta.Position,'ARPPoly')
                            % INCA UVects are DERIVED from closest approach
                            % position/velocity, not center of aperture
                            ca_pos = [polyval(output_meta.Position.ARPPoly.X(end:-1:1), ...
                                output_meta.RMA.INCA.TimeCAPoly(1)) ...
                                polyval(output_meta.Position.ARPPoly.Y(end:-1:1), ...
                                output_meta.RMA.INCA.TimeCAPoly(1)) ...
                                polyval(output_meta.Position.ARPPoly.Z(end:-1:1), ...
                                output_meta.RMA.INCA.TimeCAPoly(1))];
                            ca_vel = [polyval(polyder(output_meta.Position.ARPPoly.X(end:-1:1)), ...
                                output_meta.RMA.INCA.TimeCAPoly(1)) ...
                                polyval(polyder(output_meta.Position.ARPPoly.Y(end:-1:1)), ...
                                output_meta.RMA.INCA.TimeCAPoly(1)) ...
                                polyval(polyder(output_meta.Position.ARPPoly.Z(end:-1:1)), ...
                                output_meta.RMA.INCA.TimeCAPoly(1))];
                            if ~isfield(output_meta.RMA.INCA,'R_CA_SCP')
                                output_meta.RMA.INCA.R_CA_SCP = norm(ca_pos-SCP);
                            end
                            if (~isfield(output_meta,'Grid') || ...
                                    ((~isfield(output_meta.Grid,'Row') || ...
                                    ~isfield(output_meta.Grid.Row,'UVectECF')) && ...
                                    (~isfield(output_meta.Grid,'Col') || ...
                                    ~isfield(output_meta.Grid.Col,'UVectECF'))))
                                uRG = (SCP - ca_pos)/norm(SCP - ca_pos); % Range unit vector
                                left = cross(ca_pos/norm(ca_pos),ca_vel/norm(ca_vel));
                                look = sign(left * uRG');
                                spn=-look*cross(uRG,ca_vel); spn=spn/norm(spn); % Slant plane unit normal
                                uAZ = cross(spn,uRG);
                                output_meta.Grid.Row.UVectECF = ...
                                    cell2struct(num2cell(uRG(:)),{'X','Y','Z'});
                                output_meta.Grid.Col.UVectECF = ...
                                    cell2struct(num2cell(uAZ(:)),{'X','Y','Z'});
                            end
                        end
                        % DERIVED: Always the case for INCA
                        if isfield(output_meta,'Grid') && ...
                                isfield(output_meta.Grid,'Col') && ...
                                ~isfield(output_meta.Grid.Col,'KCtr')
                            output_meta.Grid.Col.KCtr = 0;
                        end
                        % DEFAULT: The frequency used for computing Doppler
                        % Centroid values is often the center transmitted
                        % frequency.
                        if SET_DEFAULT_VALUES && isfield(output_meta,'RadarCollection') && ...
                                isfield(output_meta.RadarCollection,'TxFrequency') && ...
                                all(isfield(output_meta.RadarCollection.TxFrequency,{'Min','Max'})) && ...
                                (~isfield(output_meta.RMA,'INCA') || ...
                                ~isfield(output_meta.RMA.INCA,'FreqZero'))
                            output_meta.RMA.INCA.FreqZero = ...
                                (output_meta.RadarCollection.TxFrequency.Min + ...
                                output_meta.RadarCollection.TxFrequency.Max)/2;
                        end
                        % Row.KCtr/FreqZero DERIVED relationship is exact
                        % (although FreqZero may be set to default above.)
                        if isfield(output_meta.RMA.INCA,'FreqZero') && ...
                                (~isfield(output_meta,'Grid') || ...
                                ~isfield(output_meta.Grid,'Row') || ...
                                ~isfield(output_meta.Grid.Row,'KCtr'))
                            output_meta.Grid.Row.KCtr = output_meta.RMA.INCA.FreqZero * 2/ SPEED_OF_LIGHT;
                        end
                    end
                end
        end
    end
end

% DERIVED: Add corners coords if they don't already exist
if (~isfield(output_meta,'GeoData')||~isfield(output_meta.GeoData,'ImageCorners')||...
        ~isfield(output_meta.GeoData.ImageCorners,'ICP'))
    try
        output_meta=add_sicd_corners(output_meta);
    end
end

% DERIVED: Add ValidData geocoords
if isfield(output_meta, 'ImageData') && isfield(output_meta.ImageData,'ValidData') && ...
        (~isfield(output_meta,'GeoData') || ~isfield(output_meta.GeoData,'ValidData'))
    try
        valid_latlons=point_slant_to_ground(...
            [[output_meta.ImageData.ValidData.Vertex.Row]; ...
            [output_meta.ImageData.ValidData.Vertex.Col]], input_meta);
        for i = 1:size(valid_latlons,2)
            output_meta.GeoData.ValidData.Vertex(i).Lat=valid_latlons(1,i);
            output_meta.GeoData.ValidData.Vertex(i).Lon=valid_latlons(2,i);
        end
    end
end
% Its difficult to imagine a scenario where GeoData.ValidData would be
% populated, but ImageData.ValidData would not, so we don't handle deriving
% the other direction.  Also, since HAE is not provided with each Vertex,
% and since the height used could have been a constant height across image
% area or from an external DEM, its not clear that there is a precise way
% to do this.

% DERIVED: Radiometric parameters RCS, sigma_0, gamma, and beta can be derived from each other
% Technically, the ratio between these terms could vary spatially across
% the image. However, this difference is usually small, so we approximate
% the relationship as a constant factor between these terms.
if isfield(output_meta, 'Radiometric') && isfield(output_meta, 'Grid') && isfield(output_meta, 'SCPCOA')
    % Calculate slant plane area
    if isfield(output_meta.Grid.Row, 'WgtFunct')
        rng_wght_f=mean(output_meta.Grid.Row.WgtFunct.^2) / ...
            (mean(output_meta.Grid.Row.WgtFunct).^2);
    else  % If no weight in metadata SICD assumes 1.0
        rng_wght_f=1.0;
    end
    if isfield(output_meta.Grid.Col, 'WgtFunct')
        az_wght_f=mean(output_meta.Grid.Col.WgtFunct.^2) / ...
            (mean(output_meta.Grid.Col.WgtFunct).^2);
    else  % If no weight in metadata SICD assumes 1.0
        az_wght_f=1.0;
    end
    area_sp = (rng_wght_f * az_wght_f) / ...
        (output_meta.Grid.Row.ImpRespBW * output_meta.Grid.Col.ImpRespBW);
    % To make the implementation shorter, first use whatever is present to
    % derive the Beta poly.
    if ~isfield(output_meta.Radiometric, 'BetaZeroSFPoly')
        if isfield(output_meta.Radiometric, 'RCSSFPoly')
            output_meta.Radiometric.BetaZeroSFPoly = ...
                output_meta.Radiometric.RCSSFPoly / area_sp;
        elseif isfield(output_meta.Radiometric, 'SigmaZeroSFPoly')
            output_meta.Radiometric.BetaZeroSFPoly = ...
                output_meta.Radiometric.SigmaZeroSFPoly / ...
                cosd(output_meta.SCPCOA.SlopeAng);
        elseif isfield(output_meta.Radiometric, 'GammaZeroSFPoly')
            output_meta.Radiometric.BetaZeroSFPoly = ...
                output_meta.Radiometric.GammaZeroSFPoly * ...
                sind(output_meta.SCPCOA.GrazeAng) / ...
                cosd(output_meta.SCPCOA.SlopeAng);
        end
    end
    % Now use the Beta poly to derive the other (if empty) fields.
    if isfield(output_meta.Radiometric, 'BetaZeroSFPoly')
        if ~isfield(output_meta.Radiometric, 'RCSSFPoly')
            output_meta.Radiometric.RCSSFPoly = ...
                output_meta.Radiometric.BetaZeroSFPoly * area_sp;
        end
        if ~isfield(output_meta.Radiometric, 'SigmaZeroSFPoly')
            output_meta.Radiometric.SigmaZeroSFPoly = ...
                output_meta.Radiometric.BetaZeroSFPoly * ...
                cosd(output_meta.SCPCOA.SlopeAng);
        end
        if ~isfield(output_meta.Radiometric, 'GammaZeroSFPoly')
            output_meta.Radiometric.GammaZeroSFPoly = ...
                output_meta.Radiometric.BetaZeroSFPoly * ...
                cosd(output_meta.SCPCOA.SlopeAng) / ...
                sind(output_meta.SCPCOA.GrazeAng);
        end
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////