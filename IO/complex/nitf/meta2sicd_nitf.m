function [ output_meta ] = meta2sicd_nitf( imgsubhdr, symmetry )
%META2SICD_NITF Converts metadata from a NITF image subheader into SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% We assume the NITF image subheader is the best source for fundamental
% fields like classification, number of rows/columns, but most other
% things we try to derive from TREs
if ~all(isspace(imgsubhdr.ISORCE)) % Can be blank
    output_meta.CollectionInfo.CollectorName=imgsubhdr.ISORCE;
end
if isempty(imgsubhdr.IID2) % The IID2 field is longer and typically a more informative image ID.
    output_meta.CollectionInfo.CoreName=imgsubhdr.IID1; % Only use IID1, if IID2 not given 
else
    output_meta.CollectionInfo.CoreName=imgsubhdr.IID2; % Use if available
end
switch imgsubhdr.ISCLAS
    case 'U'
        output_meta.CollectionInfo.Classification='UNCLASSIFIED';
    case 'S'
        output_meta.CollectionInfo.Classification='S';
    case 'T'
        output_meta.CollectionInfo.Classification='TS';
    otherwise
        output_meta.CollectionInfo.Classification='';
end
if ~isempty(imgsubhdr.ISCTLH)
    output_meta.CollectionInfo.Classification = ...
        [ output_meta.CollectionInfo.Classification '//' imgsubhdr.ISCTLH ];
end
if ~isempty(imgsubhdr.ISCODE)
    output_meta.CollectionInfo.Classification = ...
        [ output_meta.CollectionInfo.Classification '//' imgsubhdr.ISCODE ];
end

if symmetry(3)
    output_meta.ImageData.NumRows=uint32(imgsubhdr.NCOLS);
    output_meta.ImageData.NumCols=uint32(imgsubhdr.NROWS);
else
    output_meta.ImageData.NumRows=uint32(imgsubhdr.NROWS);
    output_meta.ImageData.NumCols=uint32(imgsubhdr.NCOLS);
end
output_meta.ImageData.FullImage=output_meta.ImageData; % Assume full image
output_meta.ImageData.FirstRow=uint32(0); output_meta.ImageData.FirstCol=uint32(0);

% There are multiple and redundant ways to store information in NITF.
% These TREs contain of few of the fields we are interested in.  Typically
% these will not all be found in the same image file.
if isfield(imgsubhdr,'CMETAA') % Most robust metadata for SAR data.  Use it first
    output_meta=setstructfields(meta2sicd_cmetaa(imgsubhdr.CMETAA),output_meta);
end
if isfield(imgsubhdr,'AIMIDA')
    output_meta=setstructfields(meta2sicd_aimida(imgsubhdr.AIMIDA),output_meta);
end
if isfield(imgsubhdr,'AIMIDB')
    output_meta=setstructfields(meta2sicd_aimidb(imgsubhdr.AIMIDB),output_meta);
end
if isfield(imgsubhdr,'ACFTA')
    output_meta=setstructfields(meta2sicd_acft(imgsubhdr.ACFTA),output_meta);
end
if isfield(imgsubhdr,'ACFTB')
    output_meta=setstructfields(meta2sicd_acft(imgsubhdr.ACFTB),output_meta);
end
if isfield(imgsubhdr,'BLOCKA')
    output_meta=setstructfields(meta2sicd_blocka(imgsubhdr.BLOCKA),output_meta);
end
if isfield(imgsubhdr,'EXPLTA')
    output_meta=setstructfields(meta2sicd_explt(imgsubhdr.EXPLTA),output_meta);
end
if isfield(imgsubhdr,'EXPLTB')
    output_meta=setstructfields(meta2sicd_explt(imgsubhdr.EXPLTB),output_meta);
end
% ACFTA samples spacing values are actually ground plane.  If EXPLTA/B
% gives us geometry info, we can correct to slant plane values (assuming
% our image is in the slant plane).
% If CMETAA is defined, don't bother.  Just use CMETAA values instead.
if ~isfield(imgsubhdr,'CMETAA')&&isfield(output_meta,'Grid')&&...
        ~(isfield(output_meta.Grid,'ImagePlane')&&strcmp(output_meta.Grid.ImagePlane,'GROUND'))
    if isfield(output_meta.Grid,'Col')&&isfield(output_meta.Grid.Col,'SS')&&...
            isfield(output_meta,'SCPCOA')&&isfield(output_meta.SCPCOA,'GrazeAng')
        output_meta.Grid.Col.SS=cosd(output_meta.SCPCOA.GrazeAng)*output_meta.Grid.Col.SS;
    end
    if isfield(output_meta.Grid,'Row')&&isfield(output_meta.Grid.Row,'SS')&&...
            isfield(output_meta,'SCPCOA')&&isfield(output_meta.SCPCOA,'TwistAng')
        output_meta.Grid.Row.SS=cosd(output_meta.SCPCOA.TwistAng)*output_meta.Grid.Row.SS;
    end
end
% Sensor position coordinates from MPDSRA are more reliable than MENSRA/B.
% However, MENSRA/B might be required for SCP, since this is sometimes 
% inaccurate or missing in MPDSRA.  So we process MENSR first and extract
% SCP, but we don't use other values from MENSR until end if not defined
% anywhere else.
if isfield(imgsubhdr,'MENSRA')
    mensr_meta=meta2sicd_mensra(imgsubhdr.MENSRA);
end
if isfield(imgsubhdr,'MENSRB')
    mensr_meta=meta2sicd_mensrb(imgsubhdr.MENSRB);
end
if exist('mensr_meta','var')&&...
        (~isfield(output_meta,'GeoData')||~isfield(output_meta.GeoData,'SCP'))
    output_meta.GeoData.SCP=mensr_meta.GeoData.SCP;
end
if isfield(imgsubhdr,'MPDSRA')
    ORP_ECF=[output_meta.GeoData.SCP.ECF.X...
        output_meta.GeoData.SCP.ECF.Y output_meta.GeoData.SCP.ECF.Z];
    output_meta=setstructfields(meta2sicd_mpdsra(imgsubhdr.MPDSRA,ORP_ECF),output_meta);
end
if exist('mensr_meta','var')
    output_meta=setstructfields(mensr_meta,output_meta);
end
% MPDSRA might populate PFA.FPN even for non-PFA data
if isfield(output_meta,'ImageFormation') && ...
        isfield(output_meta.ImageFormation,'ImageFormAlgo') && ...
        ~strcmp(output_meta.ImageFormation.ImageFormAlgo,'PFA') && ...
        isfield(output_meta,'PFA')
    output_meta = rmfield(output_meta,'PFA');
end

% The corner coordinates from the image subheader are generally not as
% reliable or as precise as in the above TREs (either CMETAA or BLOCKA),
% but if they are not available elsewhere, we fill them out here.
if isfield(imgsubhdr,'IGEOLO') && (~isfield(output_meta,'GeoData') || ...
        ~isfield(output_meta.GeoData,'ImageCorners') || ...
        ~isfield(output_meta.GeoData.ImageCorners,'ICP'))
    % Order of corners probably not right if symmetry~=[0 0 0]
    switch imgsubhdr.ICORDS
        case 'D'
            output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=str2double(imgsubhdr.IGEOLO(1:7));
            output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=str2double(imgsubhdr.IGEOLO(8:15));
            output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=str2double(imgsubhdr.IGEOLO(16:22));
            output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=str2double(imgsubhdr.IGEOLO(23:30));
            output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=str2double(imgsubhdr.IGEOLO(31:37));
            output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=str2double(imgsubhdr.IGEOLO(38:45));
            output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=str2double(imgsubhdr.IGEOLO(46:52));
            output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=str2double(imgsubhdr.IGEOLO(53:60));
        case 'G'
            output_meta.GeoData.ImageCorners.ICP.FRFC.Lat=latlonnum(imgsubhdr.IGEOLO(1:7));
            output_meta.GeoData.ImageCorners.ICP.FRFC.Lon=latlonnum(imgsubhdr.IGEOLO(8:15));
            output_meta.GeoData.ImageCorners.ICP.FRLC.Lat=latlonnum(imgsubhdr.IGEOLO(16:22));
            output_meta.GeoData.ImageCorners.ICP.FRLC.Lon=latlonnum(imgsubhdr.IGEOLO(23:30));
            output_meta.GeoData.ImageCorners.ICP.LRLC.Lat=latlonnum(imgsubhdr.IGEOLO(31:37));
            output_meta.GeoData.ImageCorners.ICP.LRLC.Lon=latlonnum(imgsubhdr.IGEOLO(38:45));
            output_meta.GeoData.ImageCorners.ICP.LRFC.Lat=latlonnum(imgsubhdr.IGEOLO(46:52));
            output_meta.GeoData.ImageCorners.ICP.LRFC.Lon=latlonnum(imgsubhdr.IGEOLO(53:60));
    end
end
% In general, the order of the corner coordinates pulled directly from NITF
% fields will not be correct, since each data provider has a different view
% of what the terms "upper-left", "row", "column", etc. mean with respect
% to energy down, view-from-above.  We make a generic attempt to reorient
% these based on the known orientation of the pixel data in the file with
% respect to the SICD orientation. However, if a source of the data does
% not follow this convention, this would have to be corrected in the
% sensor-specific NITF adjustments. Another (perhaps better) alternative is
% to rederive the corner coordinates using derived_sicd_fields.m (assuming
% sufficient metadata to support the SICD sensor model.)
if isfield(output_meta,'GeoData') && isfield(output_meta.GeoData,'ImageCorners') && ...
        isfield(output_meta.GeoData.ImageCorners,'ICP')
    if isequal(symmetry, [ 1 0 1 ]) % 90 degrees CCW
        temp = output_meta.GeoData.ImageCorners.ICP.FRFC;
        output_meta.GeoData.ImageCorners.ICP.FRFC = output_meta.GeoData.ImageCorners.ICP.FRLC;
        output_meta.GeoData.ImageCorners.ICP.FRLC = output_meta.GeoData.ImageCorners.ICP.LRLC;
        output_meta.GeoData.ImageCorners.ICP.LRLC = output_meta.GeoData.ImageCorners.ICP.LRFC;
        output_meta.GeoData.ImageCorners.ICP.LRFC = temp;
    elseif isequal(symmetry, [ 0 1 1 ]) % 90 degrees  CW
        temp = output_meta.GeoData.ImageCorners.ICP.FRFC;
        output_meta.GeoData.ImageCorners.ICP.FRFC = output_meta.GeoData.ImageCorners.ICP.LRFC;
        output_meta.GeoData.ImageCorners.ICP.LRFC = output_meta.GeoData.ImageCorners.ICP.LRLC;
        output_meta.GeoData.ImageCorners.ICP.LRLC = output_meta.GeoData.ImageCorners.ICP.FRLC;
        output_meta.GeoData.ImageCorners.ICP.FRLC = temp;
    elseif isequal(symmetry, [ 1 1 0 ]) % 180 degrees
        temp = output_meta.GeoData.ImageCorners.ICP.FRFC;
        output_meta.GeoData.ImageCorners.ICP.FRFC = output_meta.GeoData.ImageCorners.ICP.LRLC;
        output_meta.GeoData.ImageCorners.ICP.LRLC = temp;
        temp = output_meta.GeoData.ImageCorners.ICP.FRLC;
        output_meta.GeoData.ImageCorners.ICP.FRLC = output_meta.GeoData.ImageCorners.ICP.LRFC;
        output_meta.GeoData.ImageCorners.ICP.LRFC = temp;
    end
end
% Like corner coordinations, CollectStart is also not as reliable as in
% TREs, but if this is the only place it exist, use it.
if ~isfield(output_meta, 'Timeline') || ~isfield(output_meta.Timeline, 'CollectStart')
    try
        output_meta.Timeline.CollectStart=datenum(imgsubhdr.IDATIM,'yyyymmddHHMMSS');
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////