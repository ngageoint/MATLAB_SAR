
function pass = example_graze_pair_filt(sicdmeta1, sicdmeta2)
%EXAMPLE_GRAZE_PAIR_FILT An example image pair filter used by 'image_pair_search'
% EXAMPLE_GRAZE_PAIR_FILT(sicdmeta1, sicdmeta) takes two SICD-compatible 
% structures and returns a boolean indicating pass/fail.  Any image/SICD
% structure pair which yields a "pass" (boolean true) will be included in
% the list returned by image_pair_search.
%
% Note: No this isn't a useful filter.  At all.  You should never use it.
% We (i.e. you) should, of course, provide a filter specific to your
% particular use.  This is meant to be a stub...
%
% Just FYI...  This filter will return images whos grazing angles are
% within 5 degrees and have lat,lon locations within 0.001.  No.  There's
% no real point to it.
%
% See Also:  example_graze_filt, write_SICD_KML_placemark, 
%            write_KML_header, write_KML_footer
%
% Author: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

  pass = false;
  
  if isfield(sicdmeta1,'SCPCOA')  &&  ...
     isfield(sicdmeta1.SCPCOA,'GrazeAng')  &&  ...
     isfield(sicdmeta1, 'GeoData')  &&  ...
     isfield(sicdmeta1.GeoData, 'SCP')  &&  ...
     isfield(sicdmeta1.GeoData.SCP, 'LLH')  &&  ...
     isfield(sicdmeta2,'SCPCOA')  &&  ...
     isfield(sicdmeta2.SCPCOA,'GrazeAng')  &&  ...
     isfield(sicdmeta2, 'GeoData')  &&  ...
     isfield(sicdmeta1.GeoData, 'SCP')  &&  ...
     isfield(sicdmeta2.GeoData.SCP, 'LLH')
 
    if abs(sicdmeta1.SCPCOA.GrazeAng-sicdmeta2.SCPCOA.GrazeAng) < 5 && ...
       abs(sicdmeta1.GeoData.SCP.LLH.Lat-sicdmeta2.GeoData.SCP.LLH.Lat) < 0.001 && ...
       abs(sicdmeta1.GeoData.SCP.LLH.Lon-sicdmeta2.GeoData.SCP.LLH.Lon) < 0.001
      pass = true;
    end
    
  end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
