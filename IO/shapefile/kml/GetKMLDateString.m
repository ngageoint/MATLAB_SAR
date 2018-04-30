function outstring = GetKMLDateString( InputTime )
%GETKMLDATESTRING Formats MATLAB time into KML Date-Time string

outstring = datestr(InputTime, 'yyyy-mm-ddTHH:MM:SS.FFFZ');

end

