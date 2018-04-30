function target = text(this,long,lat,alt,txt,varargin)
%KML.TEXT(long,lat,alt,txt) Writes the text given by txt at the coordinates 
%  given by long, lat and alt. To write in more than one coordinate, pass an
%  array of coordinates, and a cell of texts, with the same number of members.
%  
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

%     target = struct('type','','id','');
    
    if ~iscell(txt)
        txt = {txt};
    end
    
    nlat = numel(lat);
    
    if ~(numel(long)==nlat && numel(alt)==nlat && numel(txt)==nlat)
        error('Invalid input sizes')
    end
         
    v = varargin;
    
    [tmpvar,hasIconURL   ]= ismemberVarargin('iconURL',v);
    [tmpvar,hasIconScale ]= ismemberVarargin('iconScale',v);
    
    if any(hasIconURL)
        v{hasIconURL+1} = 'none';
    else
        v{end+1} = 'iconURL';
        v{end+1} = 'none';
    end
    
    if any(hasIconScale)
        v{hasIconScale+1} = 0;
    else
        v{end+1} = 'IconScale';
        v{end+1} = 0;
    end
    
    for i = 1:numel(lat)
        target(i) = this.point(long(i),lat(i),alt(i),'name',txt{i},v{:});
    end
    
end

function [tf,loc] = ismemberVarargin(a,varin)
    if mod(numel(varin),2)==1
        error('Invalid number of named arguments.')
    end
    if ~iscell(a)
        a = {a};
    end
    tf  = false(size(a));
    loc = zeros(size(a));
    for i = 1:numel(a)
        for j = 1:2:numel(varin)
            if ischar(varin{j}) && strcmp(varin{j},a{i})
               tf(i) = true;
               loc(i) = j;
            end
        end
    end
end