function fftfile( compleximagefile, outfile, varargin )
%FFTFILE Computes the FFT of a file even if its too large to fit into
%memory
%
%    fftfile(compleximagefile, outfile, 'PropertyName',PropertyValue,...)
%
%       Property name     Description
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
%       sgn               Sign of the FFT to use for the transformation--
%                            -1 (default) or +1
%       dims              dimensions along which to compute FFT (options:
%                            1, 2, or [1 2] (2-D FFT), default = [1 2])
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('dims',[1 2],@(x) isequal(x,1)||isequal(x,2)||isequal(x,[1 2]));
p.addParamValue('framenumber',1);
p.addParamValue('sgn',[]);
p.FunctionName = 'fftfile';
p.parse(varargin{:});
dims=p.Results.dims;

if(length(dims)>1),
    firstoutfilename=tempname;
else
    firstoutfilename=outfile;
end
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
% The only metadata we need is the FFT sign.
% A negative FFT sign means that a forward FFT is required to transform the
% data into the spatial frequency domain.  A positive sign means an inverse
% FFT is required.  Howerver, because of the way MATLAB usually displays
% images (lowest y index at the top and increasing downward) and the way
% polar angle is typically defined (increasing counterclockwise), we flip
% the direction of the FFT to flip the direction the data is displayed.  We
% want highest RF frequency at the top of the display with columns index
% moving in the direction of a clockwise polar angle.
if ~isempty(p.Results.sgn)
    if p.Results.sgn>0
        fft_fun = {@fft, @fft};
    else % Default is -1 (fft), if not specified
        fft_fun = {@ifft, @ifft};
    end
else
    if isfield(metadata, 'Grid') && isfield(metadata.Grid,'Col') && ...
            isfield(metadata.Grid.Col,'Sgn') && metadata.Grid.Col.Sgn(1) == 1
        fft_fun{1} = @fft;
    else
        fft_fun{1} = @ifft; % Default is -1 (fft), if not specified
    end
    if isfield(metadata, 'Grid') && isfield(metadata.Grid,'Row') && ...
            isfield(metadata.Grid.Row,'Sgn') && metadata.Grid.Row.Sgn(1) == 1
        fft_fun{2} = @fft;
    else
        fft_fun{2} = @ifft; % Default is -1 (fft), if not specified
    end
end

% Process separably to file
process_by_lines( compleximagefile, firstoutfilename, ...
    @(x) fftshift(fft_fun{dims(1)}(x,[],dims(1)),dims(1)), ...
    'function_name', ['FFT dimension ' num2str(dims(1))], ...
    varargin{:}, 'dim', dims(1) );
if(length(dims)>1)
    process_by_lines( firstoutfilename, outfile, ...
        @(x) fftshift(fft_fun{dims(end)}(x,[],dims(end)),dims(end)), ...
        'function_name', ['FFT dimension ' num2str(dims(end))], ...
        'azlimits', 'full', 'rnglimits', 'full', 'dim', dims(end));
    delete(firstoutfilename);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////