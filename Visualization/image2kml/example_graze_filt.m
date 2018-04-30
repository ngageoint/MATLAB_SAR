
function pass = example_graze_filt(sicdmeta)
%EXAMPLE_GRAZE_FILT An example image filter used by 'image_search'
% EXAMPLE_GRAZE_FILT(sicdmeta) takes  a SICD-compatible structure and
% returns a boolean indicating pass/fail.  Any image/SICD structure which
% yields a "pass" (boolean true) will be included in the list returned by
% image_search.
%
% Note: No this isn't a useful filter.  At all.  You should never use it.
% We (i.e. you) should, of course, provide a filter specific to your
% particular use.  This is meant to be a stub...
%
% Just FYI...  This filter will return images whos grazing angle is between
% 30 and 40 degrees (not inclusive).  No.  There's no real point to it.
%
% Author: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

  pass = false;
  
  if isfield(sicdmeta,'SCPCOA')  &&  ...
     isfield(sicdmeta.SCPCOA,'GrazeAng')  &&  ...
     sicdmeta.SCPCOA.GrazeAng > 30  &&  ...
     sicdmeta.SCPCOA.GrazeAng < 40
    pass = true;
  end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
