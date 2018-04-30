function [ output_meta ] = meta2sicd_sio( sio_hdr, user_data )
%META2SICD_SIO Converts SIO header into SICD-style structure
%
% Really not much to do here since SIO header is so minimal.  Just fill in
% image size and datatype.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if nargin<2
    user_data=struct();
end

% Check to see if we have a SICD-in-SIO type file.  If there's a SICDMETA
% field in the user_data structure we'll assume that has all the valid
% metadata for this file and we'll use that.  Otherwise we'll just return a
% small 'stub' SICD metadata structure (not much information in the SIO
% header).
if (isfield(user_data,'SICDMETA'))
  output_meta = sicdxml2struct( xmlread( java.io.StringBufferInputStream(user_data.SICDMETA) ) );
else
  output_meta.ImageData.NumRows=sio_hdr(2);
  output_meta.ImageData.NumCols=sio_hdr(3);
  output_meta.ImageData.FullImage=output_meta.ImageData; % Assume full image, but no way to know for sure
  output_meta.ImageData.FirstRow=uint32(0); output_meta.ImageData.FirstCol=uint32(0);
  if (sio_hdr(4)==13)&&(sio_hdr(5)==8), output_meta.ImageData.PixelType='RE32F_IM32F'; end;
  output_meta.ImageData.SCPPixel.Row=fix(output_meta.ImageData.NumRows/2); % Not usually given explicitly in accompanying metadata, so just guess center
  output_meta.ImageData.SCPPixel.Col=fix(output_meta.ImageData.NumCols/2);
end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////