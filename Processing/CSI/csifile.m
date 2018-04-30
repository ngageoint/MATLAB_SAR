function csifile(compleximagefile, outfile, varargin) 
% CSIFILE Color Subaperture Image
%    csifile(compleximagefile, outfile, 'PropertyName', PropertyValue, ...)
%
% Displays subaperture information as color on full resolution data.
%
%       Property name     Description
%       dim               dimension over which to split subaperture
%                            (default = 1)
%       fill              fill factor (default = 1)
%       platformdir       platform direction, 'right' (default) or 'left'.
%                            Assumption is that 2nd dimension is increasing
%                            range.
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
%       rset             create an rset for output file (default = false)
%
% Output is stored in one color channel per file, where the first channel
% has the filename OUTFILE, and each consecutive channel has a number
% appended onto the filename.  OUTFILE = red, OUTFILE1 = green, OUTFILE2 =
% blue.
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Make editable copy of input parameters
input_params=varargin;
p = inputParser; % Extract parameter-value pairs
p.KeepUnmatched=true;
p.addParamValue('dim',1);
p.addParamValue('framenumber',1);
p.addParamValue('fill',1);
p.addParamValue('platformdir',1);
p.addParamValue('rset',false);
p.FunctionName = 'csifile';
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
        warning('CSIFILE:MISSING_FILL','No fill factor given in parameters or found in image metadata.  Assuming Nyquist-sampled data.');
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
        warning('CSIFILE:MISSING_PLATFORMDIR','No platform direction given in parameters or found in image metadata.  Assuming right looking.');
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

% Compute CSI
func_hand=@(x) mat2cell(csimem(x, input_params{:}),size(x,1),size(x,2),[1 1 1]);
process_by_lines( compleximagefile, outfile, func_hand, ...
    'function_name', 'CSI', input_params{:} );
if apply_normalization, delete(compleximagefile); end; % Remove temporary file

if exist(outfile,'file')
    % Create wrapper file for all three color bands
    fid=fopen([outfile '.mbw'],'w');
    fprintf(fid,'MULTI-BAND WRAPPER\n');
    fprintf(fid,'TYPE=RGB\n');
    fprintf(fid,'R=%s\n', outfile);
    fprintf(fid,'G=%s002\n', outfile);
    fprintf(fid,'B=%s003', outfile);
    fclose(fid);
    if p.Results.rset
        myImageAdapter = ComplexSarRemapAdapter([outfile '.mbw']);
        rsetwrite(myImageAdapter,[outfile '.rset']);
        myImageAdapter.close();
        delete([outfile '.mbw']);
        delete(outfile);
        delete([outfile '002']);
        delete([outfile '003']);
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////