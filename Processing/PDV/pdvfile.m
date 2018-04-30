function pdvfile( compleximagefile, outfile, varargin )
%PDVFILE Phase Derivative Value
%    pdvfile(compleximagefile, outfile, 'PropertyName',PropertyValue,...)
%    calculates the phase derivative value of compleximagefile using the
%    properties specified. This version of PDV does NOT require that the
%    complete data fit into memory.  It processes from any format handled
%    by OPEN_READER and outputs to files in SIO format.
%
%       Property name     Description
%       deltax            pixel shift (default = 0.25)
%       filtersize        size of smoothing filter (default = [5 5])
%       filtertype        type of filter, 'mean' (default) or 'median'
%       dim               dimension over which to calculate phase
%                            derivative (default = 1)
%       azlimits          min and max samples in azimuth over which to
%                            compute (default = query user with GUI).  Use
%                            'full' to specify entire range.  If nothing is
%                            selected, a MATLAB GUI will be used to select
%                            an area.
%       rnglimits         min and max lines in range over which to compute
%                            (default = query user with GUI).  Use 'full'
%                            to specify entire range.  If nothing is
%                            selected, a MATLAB GUI will be used to select
%                            an area.
%       framenumber      frame to process if a multi-image file
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Get image metadata
func_args=parsepdvinputs(varargin{:});
input_params=varargin;

% Read image metadata
fr=open_reader(compleximagefile);
if iscell(fr)
    metadata=fr{func_args.framenumber}.get_meta();
    for i=1:length(fr)
        fr{i}.close();
    end
else
    metadata=fr.get_meta();
    fr.close();
end

try
    % Normalize data if necessary.
    apply_normalization = ~is_normalized_sicd(metadata,func_args.dim);
    if apply_normalization
        tmp=tempname(); % Stored deskewed data in temporary file to process
        normalize_complex_file(compleximagefile,'output_filename',tmp,varargin{:});
        compleximagefile=tmp;
        input_params=[input_params {'azlimits','full','rnglimits','full'}];
    end
catch
    apply_normalization = 0;
end

% Actual PDV
func_hand=@(x) pdvmem(x, func_args);
process_by_lines( compleximagefile, outfile, func_hand,...
    'block_overlap', floor(func_args.filtersize.*[1 1]/2),...
    'function_name', 'PDV', input_params{:} );
if apply_normalization, delete(compleximagefile); end; % Remove temporary file

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////