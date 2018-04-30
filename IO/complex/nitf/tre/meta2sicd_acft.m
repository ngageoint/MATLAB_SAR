function [ output_meta ] = meta2sicd_acft( acft_struct )
%META2SICD_ACFT Converts metadata stored in the ACFTA/B NITF TRE into a
% SICD metadata format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

M2F=3.28083989501312; % Set meters-to-feet conversion constant

if ~isempty(acft_struct.SENSOR_ID)
    output_meta.CollectionInfo.CollectorName=acft_struct.SENSOR_ID;
end

% This NITF value is actually the ground plane spacing.  Should actually be
% multiplied by cos of graze angle for SICD slant plane spacing.
output_meta.Grid.Col.SS=acft_struct.ROW_SPACING;
% This NITF value is actually the ground plane spacing.  Should actually be
% multiplied by cos of tilt angle for SICD slant plane spacing.
output_meta.Grid.Row.SS=acft_struct.COL_SPACING;

if ~isfield(acft_struct,'ROW_SPACING_UNITS')||acft_struct.ROW_SPACING_UNITS=='f'
    output_meta.Grid.Col.SS=output_meta.Grid.Col.SS/M2F;
end
if ~isfield(acft_struct,'COL_SPACING_UNITS')||acft_struct.COL_SPACING_UNITS=='f'
    output_meta.Grid.Row.SS=output_meta.Grid.Row.SS/M2F;
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////