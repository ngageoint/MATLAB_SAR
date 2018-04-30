function [ readerobj ] = open_sio_reader( filename )
%OPEN_SIO_READER Intiates a reader object for SIO files
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Setup reader type
[ihdr,endian,data_offset,user_data] = read_sio_meta( filename );
meta=meta2sicd_sio(ihdr,user_data); % Convert image size info to SICD struct
meta.native.sio=user_data;

%% Check for CASPR metadata file and read
caspr_filename=locate_caspr(filename);
symmetry=[0 0 0]; % Energy from top
if ~isempty(caspr_filename)
    native_metadata=read_caspr_meta(caspr_filename);
    if isfield(native_metadata,'ImageParameters')&&...
            isfield(native_metadata.ImageParameters,'imageilluminationdirectiontopleftbottomright')
        if strcmpi(native_metadata.ImageParameters.imageilluminationdirectiontopleftbottomright,'left')
            symmetry=[0 1 1];
            meta=meta2sicd_sio(ihdr([1 3 2 4 5])); % Reorient size info
        elseif ~strcmpi(native_metadata.ImageParameters.imageilluminationdirectiontopleftbottomright,'top')
            error('OPEN_SIO_READER:UHANDLED_ILLUMINATION_DIRECTION','Unhandled illumination direction.');
        end
    end
    meta=setstructfields(meta2sicd_caspr(native_metadata),meta);
    meta.native.caspr=native_metadata;
    
    meta = add_sicd_corners(meta); % Corners coords can only be computed with both SIO and CASPR metadata
end

%% Build object
datasize=ihdr([3 2]);
[datatype,complexbool,freadtype]=siotype2matlab(ihdr(4),ihdr(5));
bands=1;
readerobj=open_generic_reader(filename, datasize, freadtype, complexbool,...
    data_offset, endian, symmetry, bands);
readerobj.get_meta=@(x) meta;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////