function subaperturefile( compleximagefile, outfile, varargin )
%SUBAPERTUREFILE Subaperture processing on a file
%    subaperturefile(compleximagefile, outfile, 'PropertyName',PropertyValue,...)
%
% Calculates the subapertures of compleximagefile using the properties
% specified. This version of subaperture does NOT require that the complete data fit
% into memory.  It processes from any format handled by OPEN_READER and
% outputs to files in SIO format.
%
%       Property name     Description
%       frames            number of frames (default = 7)
%       apfraction        fraction of aperture for each subaperture
%                            (default = .25)
%       method            'normal' (default), 'fullpixel', or 'minimal'
%       platformdir       platform direction, 'right' (default) or 'left'
%       dim               dimension over which to split subaperture
%                            (default = 1)
%       fill              fill factor
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
% Output frames are stored in one frame per file, where the first frame has
% the filename OUTFILE, and each consecutive frame has the frame number
% appended onto the filename.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Make editable copy of input parameters
input_params=varargin;
p = inputParser; % Extract parameter-value pair required for difile
p.KeepUnmatched=true;
p.addParamValue('dim',1);
p.addParamValue('framenumber',1);
p.addParamValue('fill',1);
p.addParamValue('platformdir',1);
p.FunctionName = 'subaperturefile';
p.parse(varargin{:});

% Get image metadata
fr=open_reader(compleximagefile);
if iscell(fr)
    metadata=fr{p.Results.framenumber}.get_meta();
    for i=1:length(fr)
        fr{i}.close();
    end
else
    metadata=fr.get_meta();
    fr.close();
end
% Get metadata to calculate fill factor
if any(strcmp(p.UsingDefaults,'fill'))
    if isfield(p.Results,'dim')&&(p.Results.dim==2)&&...
            isfield(metadata,'Grid')&&isfield(metadata.Grid,'Row')&&...
            all(isfield(metadata.Grid.Row,{'SS','ImpRespBW'}))
        input_params=[{'fill' 1/(metadata.Grid.Row.SS*metadata.Grid.Row.ImpRespBW)} input_params]; % Rng Zeropad
    elseif isfield(metadata,'Grid')&&isfield(metadata.Grid,'Col')&&...
            all(isfield(metadata.Grid.Col,{'SS','ImpRespBW'}))
        input_params=[{'fill' 1/(metadata.Grid.Col.SS*metadata.Grid.Col.ImpRespBW)} input_params]; % Az Zeropad
    else
        warning('SUBAPERTUREFILE:MISSING_FILL','No fill factor given in parameters or found in image metadata.  Assuming Nyquist-sampled data.');
    end
end
% Platform direction
if any(strcmp(p.UsingDefaults,'platformdir'))
    if isfield(metadata,'SCPCOA')&&isfield(metadata.SCPCOA,'SideOfTrack')
        if metadata.SCPCOA.SideOfTrack=='R'
            input_params=[{'platformdir' 'right'} input_params];
        elseif metadata.SCPCOA.SideOfTrack=='L'
            input_params=[{'platformdir' 'left'} input_params];
        end
    else % Check for platformdir in varargin; otherwise return warning
        warning('SUBAPERTUREFILE:MISSING_PLATFORMDIR','No platform direction given in parameters or found in image metadata.  Assuming right looking.');
    end
end

% Normalize data if necessary.
apply_normalization = ~is_normalized_sicd(metadata,p.Results.dim);
if apply_normalization
    tmp=tempname(); % Stored deskewed data in temporary file to process
    normalize_complex_file(compleximagefile,'output_filename',tmp,varargin{:});
    compleximagefile=tmp;
    input_params=[input_params {'azlimits','full','rnglimits','full'}];
end

% Parse input parameters
func_args=parsesubapertureinputs(input_params{:});

% Compute subapertures
func_hand=@(x) subaperturemem(x, func_args);
process_by_lines( compleximagefile, outfile, func_hand, ...
    'function_name', 'subapertures', input_params{:} );
if apply_normalization, delete(compleximagefile); end; % Remove temporary file

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////