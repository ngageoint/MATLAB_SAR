classdef kml < handle
%KML(name) Create a KML toolbox object, that allows you to plot from MATLAB to Google Earth
%  This class-based toolbox allows you to create many different plots in Google Earth, by automatically creating the required xml-based kml files without user interaction.
%  
%  With it, you can create: 
%  - line plots, scatter plots 
%  - 2D and 3D contours 
%  - 2D and 3D polygons 
%  - quiver plots 
%  - write text in a given point 
%  - place 3D models 
%  - overlay images 
%  - transfer more complex figures as images 
%  - folders, subfolders,... to aggregate similar plots 
%  
%  For some usage examples, check the code at .\test\RunTests.m
%
%  If you enjoy it, just drop me an email at kml@rafael.aero saying for what you're using it :).
%    
%   Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)
%   $Revision: 2.3 $  $Date: 2012/09/05 08:00:00 $

    properties %(Access = private, Transient)
        xml
        doc
        parent = []
        filename = ''
        name = ''
        unit = 'deg'
        zip = true
        includeFiles;
    end
    
    properties (Access = private)
        target = 'earth'
    end
    
    methods
        function this = kml(name,parent)
            if nargin <1
                name = 'Unnamed KML';
            end
            this.name = name;
            if nargin <2
                parent = [];
            else
                if ~isa(parent,'kml')
                    error('Invalid parent, pass only a kml object or a kml folder')
                end
            end
            
            this.parent = parent;
            
            if isempty(parent)
                this.xml = com.mathworks.xml.XMLUtils.createDocument('kml');
                kmlNode = this.xml.getDocumentElement;
                kmlNode.setAttribute('xmlns','http://www.opengis.net/kml/2.2');
                kmlNode.setAttribute('xmlns:gx','http://www.google.com/kml/ext/2.2');
                kmlNode.setAttribute('xmlns:kml','http://www.opengis.net/kml/2.2');
                kmlNode.setAttribute('xmlns:atom','http://www.w3.org/2005/Atom');

                this.doc = this.xml.createElement('Document');

                this.doc.appendChild(this.textNode('name',name));
            
                kmlNode.appendChild(this.doc);
            else
                this.doc = parent.xml.createElement('Folder');
                this.xml = parent.xml;
                
                this.doc.appendChild(this.textNode('name',name));
                this.doc.appendChild(this.textNode('id',  name));
                
                parent.doc.appendChild(this.doc);
            end
        end
        
        function clear(this)
            if isempty(this.parent)
                %not a folder, just create a new kml class and reuse it's xml nodes
                a = kml(this.name);
                this.doc = a.doc;
                this.xml = a.xml;
                this.includeFiles = '';
                a.delete;
            else
                this.parent.doc.removeChild(this.doc);
                a = kml(this.name,this.parent);
                this.doc = a.doc;
                this.xml = a.xml;
                a.delete;
            end
        end
        
        function addIncludeFile(this,varargin)
            K = this;
            while ~isempty(K.parent)
                K = K.parent;
            end
            newFiles = varargin;
            for i = 1:numel(newFiles)
                K.includeFiles{end+1} = newFiles{i};
            end
        end
        
        
        function save(this,filename, compact)
            if ~isempty(this.parent)
                warning('KML:incompletefile','You are saving only a folder of the KML file, not the full KML, this will not load on any viewer');
            end
            
            if nargin < 2
                if isempty(this.filename)
                    if strcmpi(this.name,'Unnamed KML')
                        [tmpvar,filename] = fileparts(tempname(cd));
                        filename = [filename '.kml'];
                        warning('KML:tempfile','Please delete file %s later, if not required...',filename);
                    else
                        filename = this.name;
                    end
                else
                    filename = this.filename;
                end
            end
            
            if nargin < 3
                compact = this.zip;
            else
                this.zip = compact;
            end
            
            [p,fn,ext] = fileparts(filename);
            filename = fullfile(p,[fn '.kml']);
            
            
            xmlwrite(filename,this.xml); 
            
            if this.zip
                [p,fn,ext] = fileparts(filename);
                kmzFilename = fullfile(p,[fn '.kmz']);
                zip(kmzFilename,vertcat({filename},this.includeFiles(:)).');
                delete(filename);
                movefile([kmzFilename '.zip'],kmzFilename);
                this.filename = kmzFilename;
            else
                this.filename = filename;
            end
        end
        
        function viewKML(this)
            tmpfile = [tempname '.kml'];
            xmlwrite(tmpfile ,this.xml);
            edit(tmpfile);
            delete(tmpfile);
        end
        
        function run(this)
            this.save;
            switch true
                case ispc
                    system(['start "kmltoolbox" "' this.filename '"']);
                case ismac
                    system(['open' ' "' this.filename '"']);
            otherwise
                disp('The KML file has been saved, open it in Google Earth');
            end
        end

        function useRadians(this)
            this.unit = 'rad';
        end

        function useDegrees(this)
            this.unit = 'deg';
        end        
        
        function setTarget(this,target)
            target = lower(target);
            assert(ismember(target,{'earth','moon','sky','mars'}),'Target must be one of the following: ''earth'' ''moon'' ''mars'' or ''sky''')
            this.target = target;
            kmlNode = this.xml.getDocumentElement;
            if ~strcmpi(this.target,'earth')
                kmlNode.setAttribute('hint',['target=' this.target]);
            end            
        end
        
        function varargout = checkUnit(this,varargin)
           varargout = varargin;
           
           parent = this;
           
           while ~isempty(parent.parent)
               parent = parent.parent;
           end
           
           if strcmpi(parent.unit,'deg')
               %do nothing
           elseif strcmpi(parent.unit,'rad')
               for i = 1:numel(varargin)
                   varargout{i} = varargin{i}.* (180/pi);
               end
           else
               error('invalid angular unit %s', parent.unit);
           end
        end
        
        function disp(this)
            display(this);
        end
        
        
        function display(this)
            if getpref('kmltoolbox','ShowDisclaimer',true)
                disp(sprintf('\n'));
                disp('     _  ____  __ _      _____         _ _             ')
                disp('    | |/ /  \/  | |    |_   _|__  ___| | |__  _____ __')
                disp('    | '' <| |\/| | |__    | |/ _ \/ _ \ | ''_ \/ _ \ \ /')
                disp('    |_|\_\_|  |_|____|   |_|\___/\___/_|_.__/\___/_\_\')
                disp(sprintf('\n'));
                disp('Copyright 2012 Rafael Fernandes de Oliveira (rafael@rafael.aero)')
                disp(sprintf('\n'));
                disp('If you enjoy it, just drop me an email at kml@rafael.aero saying for what you''re using it :-)')
                disp(sprintf('\n'));
                disp('To install the toolbox, and suppress this message, run kml.install')
            end
        end
        
    end
    
    methods %(Access = private)
        function node = textNode(this,Name,Text)
            node = this.xml.createElement(Name);
            node.appendChild(this.xml.createTextNode(Text));
        end
    end
    
    methods (Static)
       iconURL = parseIconURL(iconURL) 
       im      = captureScreen(view)
       
       function colorHex  = color2kmlHex(color)
           if numel(color) ==3
               color(4) = 1;
           end
           color = min(max(floor(color*255),0),255);
           [r,g,b,a] = deal(color(1),color(2),color(3),color(4));
           [rhex, ghex, bhex, ahex ]= deal(dec2hex(r),dec2hex(g),dec2hex(b),dec2hex(a));
           if length(rhex)==1,rhex=['0' rhex];end
           if length(ghex)==1,ghex=['0' ghex];end
           if length(bhex)==1,bhex=['0' bhex];end
           if length(ahex)==1,ahex=['0' ahex];end
           
           colorHex = [ahex bhex ghex rhex];
       end
       
       function ID = getTempID(rootID)
           [tmp,tmpID] = fileparts(tempname);
            ID = [rootID tmpID];
       end

        
       function dist = ll2dist(lon1, lat1,lon2, lat2)
           %Input in radians - Calculate distance using Haversine formula
           R = 6378137;
           dLat = (lat2-lat1);
           dLon = (lon2-lon1);
           a = sqr(sin(dLat/2)) + cos(lat1) .* cos(lat2) .* sqr(sin(dLon/2));
           c = 2 * atan(sqrt(a)./sqrt(1-a));
           dist = R * c;
       end
       
       function install()
           clc
           nl = sprintf('\n');
           disp(nl);
           disp('     _  ____  __ _      _____         _ _             ')
           disp('    | |/ /  \/  | |    |_   _|__  ___| | |__  _____ __')
           disp('    | '' <| |\/| | |__    | |/ _ \/ _ \ | ''_ \/ _ \ \ /')
           disp('    |_|\_\_|  |_|____|   |_|\___/\___/_|_.__/\___/_\_\')
           disp(nl);
           disp('This will add the KML toolbox to your MATLAB path, at the current location');
           disp('If you prefer to have it located somewhere else, first copy the whole folder');
           disp('to the place you want, and then run kml.install from there!');
           disp(nl);
           p = fileparts(fileparts(mfilename('fullpath')));
           disp(sprintf('Toolbox Path: %s',p));
           disp(nl);
           s = input('Continue [Y]/N?   ','s');
           if strcmpi(s,'n')
               disp('Instalation aborted!')
               return
           end
           disp(nl);
           disp('Adding toolbox folder to the MATLAB path...')
           disp(nl);
           disp('(In some systems, this might require administrative rights, sorry for that)')
           disp(nl);
           
           addpath(p);
           addpath(fullfile(p,'html'));
           r = savepath;
           
           if r == 1
               disp('Error saving your path. Please add the toolbox folder manually to your MATLAB path')
           end
           
           disp('Installing the toolbox help...')

           warning('off','MATLAB:doc:DocNotInstalled')
           builddocsearchdb(fullfile(p,'html'))
           warning('on','MATLAB:doc:DocNotInstalled')
           
           setpref('kmltoolbox','ShowDisclaimer',false);
           
           clc
           disp(nl);
           disp('     _  ____  __ _      _____         _ _             ')
           disp('    | |/ /  \/  | |    |_   _|__  ___| | |__  _____ __')
           disp('    | '' <| |\/| | |__    | |/ _ \/ _ \ | ''_ \/ _ \ \ /')
           disp('    |_|\_\_|  |_|____|   |_|\___/\___/_|_.__/\___/_\_\')
           disp(nl);           
           disp(nl);
           disp('    Toolbox installed successfully!')
           disp(nl);
           disp('To take a look in the toolbox contents, run: >> kmldoc');
           disp('To get help for the toolbox functions, run:  >> kmldoc kml.FunctionName');
           disp(nl);
           
       end
    end
end