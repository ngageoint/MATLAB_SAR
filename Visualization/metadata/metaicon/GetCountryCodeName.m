function OutString = GetCountryCodeName(in,OutFormat)
%GETCOUNTRYCODENAME Returns CC or Name based on input Lat/Lon or CC

%country lookup for lat/lon input is done via TZ world shapefile

%valid output types are: 'Digraph','Trigraph' and 'Name'

OutString = [];
Country = 0;

%determine input format
if isnumeric(in) && length(in) == 2
   %lat/lon
   InType = 'LatLon';
elseif ischar(in) && length(in) == 2
   %digraph
   InType = 'Digraph';
elseif ischar(in) && length(in) == 3
   %trigraph
   InType = 'Trigraph';
else
   return; 
end

%load csv file that has TZID,Digraph,Trigraph and CountryName
DataDirstr = which('DataDir');    
pathstr = fileparts(which('DataDir'));
fid = fopen([pathstr filesep 'TZ_CC_NAMES.csv']);
tline = fgetl(fid);
tline = fgetl(fid);
count = 0;
while ischar(tline)
    count = count + 1;
    st = strsplit(tline,',');
    TZIDs{count} = st{1};
    Digraphs{count} = st{2};
    Trigraphs{count} = st{3};
    Names{count} = st{4};
    tline = fgetl(fid);
end
fclose(fid);

%now determine country based on intype
switch InType
    case 'LatLon'
        %load shapefile to see if lat/lon is inside a polygon
        tzname = [pathstr filesep 'TimeShapes' filesep 'tz_world_mp.shp'];
        s = shaperead(tzname);
        Lat = in(1);
        Lon = in(2);
        TZID = [];
        for i=1:length(s)
            if inpolygon(Lon,Lat,s(i).X,s(i).Y)
                TZID = s(i).TZID;
                break;
            end
        end
        for i=1:count
            %try to match digraph
            if strcmpi(TZID,TZIDs{i})
                Country = i;
                break;
            end
        end
    case 'Digraph'
        for i=1:count
            %try to match digraph
            if strcmpi(in,Digraphs{i})
                Country = i;
                break;
            end
        end
    case 'Trigraph'
        for i=1:count
            %try to match digraph
            if strcmpi(in,Trigraphs{i})
                Country = i;
                break;
            end
        end
end

if Country == 0
    %no country found
    return;
end

%determine output based on country index and output format
switch OutFormat
    case 'Digraph'
        OutString = Digraphs{Country};
    case 'Trigraph'
        OutString = Trigraphs{Country};
    case 'Name'        
        OutString = Names{Country};
end


end

