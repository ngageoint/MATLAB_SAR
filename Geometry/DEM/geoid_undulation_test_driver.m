
function geoid_undulation_test_driver
% GEOID_UNDULATION_TEST_DRIVER is a test driver for the geoid_undulation
% function.  It compares computed Geoid undulation values to "reference"
% values as found in the "output.dat" file as found at
%
% http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
%
% See also: geoid_undulation
%
% Written by: Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

    filename = 'Und_min2.5x2.5_egm2008_isw=82_WGS84_TideFree_SE';

    % These values are the expected results.  The first 6 are from the
    % documentation "output.dat" file while the last 2 are "interesting
    % points" whose results have been computed via the NGA web site.  Note,
    % we don't expect to get _exactly_ the same answer for the last two
    % since the NGA site probably uses the 1' grid rather than the 2.5'
    % grid.  The results should be close though.
    correct = [ 37.000000  241.000000    -26.151;...
                37.000000 -119.000000    -26.151;...
                36.000000  242.983333    -29.171;...
                90.000000    0.000000     14.899;...
               -90.000000  359.983333    -30.150;...
               -90.000000    0.000000    -30.150;...
                42.244169 -83.5066427    -34.89;...
                33.276390 -112.282220    -30.86];
    fprintf(1,'\n');
    fprintf(1,'\n');
    fprintf(1, '---------------------------------------------------------------\n');
    fprintf(1, '    lat          lon       correct   computed   error\n');
    fprintf(1, '---------------------------------------------------------------\n');
    for i=1:length(correct)
        lat = correct(i,1);
        lon = correct(i,2);
        H   = geoid_undulation(lat, lon, 'useHiRes', true);
        
        % special handling...
        if (i==3)
            spec = '?';
        elseif (i==7  ||  i==8)
            spec = '*';
        else
            spec = ' ';
        end
            
        fprintf(1, '%10.6f %12.6f %9.3f  %9.3f   %6.3f %s\n', ...
               lat, lon, correct(i,3), H, H-correct(i,3), spec);
    end
    fprintf(1,'\n');
    fprintf(1,'? Not sure why this is "wrong".  Typo? Rounding?\n');
    fprintf(1,'* We don''t expect these to be perfect, only "close".\n');
    fprintf(1,'\n');
end
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
