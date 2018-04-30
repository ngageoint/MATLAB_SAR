function savetovariable = chipfile( compleximagefile, outfile, varargin )
%CHIPFILE Saves an AOI from a complex image to file or variable
%
%    chipfile(compleximagefile, outfile, 'PropertyName',PropertyValue,...)
%
%       Property name     Description
%       azlimits          Min and max samples in azimuth over which to
%                            compute (default = query user with GUI).  Use
%                            'full' to specify entire range.  If nothing is
%                            selected, a MATLAB GUI will be used to select
%                            an area.
%       rnglimits         Min and max lines in range over which to compute
%                            (default = query user with GUI).  Use 'full'
%                            to specify entire range.  If nothing is
%                            selected, a MATLAB GUI will be used to select
%                            an area.
%       framenumber       Frame to process if a multi-image file.  Default
%                            is 1.
%       file_format       Format of output file.  'SIO', 'SICD', or 'NRL'.
%                            Default is SICD.
%       block_size        Determines how many bytes are in memory at once.
%                            Set to Inf to read whole image at once.
%                            Default is about 50 Mbytes.
% 
% Optionally one can also use this function to chip a file to a variable in
% memory:
%    savetovariable = chipfile(compleximagefile);
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Select input image through GUI if not provided (and save path)
if (nargin<1)||isempty(compleximagefile)
    % Recall last interactively selected path used
    if ispref('matlab_sar_toolbox','last_used_directory')
        pathname = getpref('matlab_sar_toolbox','last_used_directory');
        if ~ischar(pathname)||~exist(pathname,'dir')
            pathname = pwd;
        end
    else
        pathname = pwd;
    end
    [compleximagefile,pathname]=uigetfile(sar_file_extensions('complex'),...
        'Select Input File',pathname);
    if(~compleximagefile), return; end; % Cancel was chosen
    setpref('matlab_sar_toolbox','last_used_directory',pathname);
    compleximagefile=fullfile(pathname,compleximagefile);
else
    pathname=fileparts(which(compleximagefile));
end

% Select output image through GUI if not provided (and save path)
default_filetype = 'SICD'; % If not chosen through file GUI
if ((nargin<2)||isempty(outfile))&&(nargout==0)
    [outfile,pathname,filterindex]=uiputfile(...
        {'*.ntf','SICD (*.ntf)'; '*.sio','SIO (*.sio)'; '*.iq8','NRL (*.iq8)'},...
        'Select output filename',pathname);
    if(~outfile), return; end; % Cancel was chosen
    exts = {'SICD','SIO','NRL'};
    default_filetype = exts{filterindex};
    outfile=fullfile(pathname,outfile);
end

% Parse input parameters
p = inputParser;
p.KeepUnmatched=true;
p.addParamValue('azlimits','default');
p.addParamValue('rnglimits','default');
p.addParamValue('framenumber',1);
p.addParamValue('file_format',default_filetype);
p.addParamValue('block_size',2^26); % in bytes
p.FunctionName = mfilename;
p.parse(varargin{:});
azlimits=p.Results.azlimits;
rnglimits=p.Results.rnglimits;
framenumber=p.Results.framenumber;
if(strcmpi(azlimits,'default')&&strcmpi(rnglimits,'default'))
    [aoi_info, framenumber]=mitm_viewer(compleximagefile,'mode','aoi',...
        'closeAfterSelect',true);
    if ~isempty(aoi_info) % User may have closed viewer without selecting AOI
        azlimits=[aoi_info(1) aoi_info(1)+aoi_info(3)-1];
        rnglimits=[aoi_info(2) aoi_info(2)+aoi_info(4)-1];
    else
        savetovariable = []; return;
    end
end
reader_object=open_reader(compleximagefile);
if iscell(reader_object)
    reader_object = reader_object{framenumber};
end
meta = reader_object.get_meta();
if(isempty(azlimits)||any(strcmpi(azlimits,{'full','default'})))
    azlimits=[1 double(meta.ImageData.NumCols)];
end
if(isempty(rnglimits)||any(strcmpi(rnglimits,{'full','default'})))
    rnglimits=[1 double(meta.ImageData.NumRows)];
end
save_to_file=exist('outfile','var')&&~isempty(outfile); % Will be saving chip to file
if save_to_file
    if ~isfield(meta.ImageData, 'PixelType')
        meta.ImageData.PixelType = 'RE32F_IM32F'; % Default type
    end
    % Update SICD coordinates to reflect chip area within larger image
    meta.ImageData.NumRows = uint32(diff(rnglimits)+1);
    meta.ImageData.NumCols = uint32(diff(azlimits)+1);
    meta.ImageData.FirstRow = uint32(meta.ImageData.FirstRow + rnglimits(1) - 1);
    meta.ImageData.FirstCol = uint32(meta.ImageData.FirstCol + azlimits(1) - 1);
    meta = add_sicd_corners(meta); % Update corners coords
    switch p.Results.file_format
        case 'SIO'
            writer_object=SIOWriter(outfile, meta);
        case 'SICD'
            writer_object=SICDWriter(outfile, meta);
        case 'NRL'
            writer_object=NRLWriter(outfile, meta);
        otherwise
            % Error
    end
end
if nargout>0, savetovariable=zeros(diff(azlimits)+1,diff(rnglimits)+1); end; % Saving to variable in memory

% Do chipping
first_row_in_block=rnglimits(1);
num_rows_in_block=max(1,floor(p.Results.block_size/((diff(azlimits)+1)*8))); % Assume 8-byte elements (largest allowed in SICD)
wb_hand=waitbar(0,'Chipping data...');
while(first_row_in_block<rnglimits(2))
    row_indices=first_row_in_block:min(rnglimits(2),first_row_in_block+num_rows_in_block-1);
    image_block = reader_object.read_chip(azlimits,row_indices([1 end]));
    if nargout>0, savetovariable(:,row_indices-rnglimits(1)+1)=image_block; end;
    if save_to_file
        writer_object.write_chip(image_block, [1 (first_row_in_block-rnglimits(1)+1)]);
    end
    first_row_in_block=row_indices(end)+1;
    waitbar(double(first_row_in_block-rnglimits(1))/double(diff(rnglimits)+1),wb_hand);
end
close(wb_hand);
reader_object.close();
% writer_object will close itself upon leaving this function

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////