function AddShape(fid,format,shape)
%ADDSHAPE Writes shape to shapefile of specified format
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Function appends shape to Google KML file, Remoteview RVL file, or Socet
% GRFX file.  Expects an open file handle to the shapefile (use
% OpenShapefile.m to open file).  Shape must be a Matlab structure that
% contains the following parameters:
%
%   ShapeType     - String - Enumerated Type of shape. Options are:
%                            'POINT','LINE','POLYGON','OVERLAY'
%   Lat           - Float  - Array of Latitude Values (deg)
%   Lon           - Float  - Array of Longitude Values (deg)
%   Color         - String - Shape Color (Red,Yellow,Cyan,Blue,Green,Magenta,White,Black or RGB Triplet)
%   LineThickness - Int    - Line Thickness (1-5) 
%   FillColor     - String - Fill Color (Red,Yellow,Cyan,Blue,Green,Magenta,White,Black)
%   Transparency  - Float  - Transparency for fill (1=Opaque,0=Transparent)
%   Name          - String - Shape Name
%
% KML file
%   Alt           - Float  - Altitude (meters above MSL).  Optional
%   Date          - String - Date Stamp.  Must be in KML Format (i.e.
%                            YYYY-MM-DD).  If empty, then ignored.
%   StartTime     - String - Start Time Stamp. Must also include a Stop
%                            Time. Format is YYYY-MM-DDTHH:MM:SSZ
%   StopTime      - String - Stop Time Stamp. Must also include a Start
%                            Time. Format is YYYY-MM-DDTHH:MM:SSZ
%   Description   - String - Shape Description (html used in kml placemark)
%   OverlayBox    - Float  - Array of Overlay Bounding Box (required for
%                            overlay type) [north,south,east,west]
%
% SocetGXP GRFX File
%   VertexFlag    - Int    - Vertex Flag (1=Display,0=Hide)
%   VertexSize    - Int    - Vertex Size (1-5)
%   VertexColor   - String - Vertex Color (Red,Yellow,Cyan,Blue,Green,Cyan,Magenta,White,Black)
%   LineFlag      - Int    - Line Flag (1=Display, 0=Hide)
%   LineColor     - String - Line Color (Red,Yellow,Cyan,Blue,Green,Cyan,Magenta,White,Black)
%
%INPUTS:
%  fid      - required: shapefile fid
%  format   - required: 'grfx','rvl' or 'kml'
%  shape    - required: structure with shape specification
%
%OUTPUTS:
% 
%VERSION:
%   1.0 
%     - Tim Cox 20090324
%     - initial version
%   1.1
%     - Tim Cox 20120514
%     - Added altitude for KML output

format = lower(format);

%if required fields are not in structure then return
if (~isfield(shape,'ShapeType'))
    error('ADDSHAPE:NOTYPE', 'Cannot add shape without shape type');
end
if (~isfield(shape,'Lat') || ~isfield(shape,'Lon')) && ~strcmp(shape.ShapeType,'OVERLAY')
    error('ADDSHAPE:NOLATLON', 'Cannot add shape without Lat/Lon Data');
end

%check for Altitude.  If not specified then set to zero
if (~isfield(shape,'Alt')) && ~strcmp(shape.ShapeType,'OVERLAY')
    shape.Alt = zeros(length(shape.Lat),1);
    AltFlag = 0;
else
    AltFlag = 1;
end

%set any fields that are not in the shape structure to default values
if (~isfield(shape,'LineThickness'))
    shape.LineThickness = 1;
end
if (~isfield(shape,'Transparency'))
    shape.Transparency = .5;
end
if (~isfield(shape,'Color'))
    shape.Color = 'Red';
end
if (~isfield(shape,'FillFlag'))
    shape.FillFlag = isfield(shape,'FillColor')&&~isempty(shape.FillColor);
end
if (~isfield(shape,'FillColor'))
    shape.FillColor = shape.Color;
end
if (~isfield(shape,'Name'))
    shape.Name = 'Shape';
end
if (~isfield(shape,'Description'))
    shape.Description = '';
end

%grfx only fields
if (strcmp(format,'grfx') == 1)
    if (~isfield(shape,'VertexFlag'))
        shape.VertexFlag = 0;
    end
    if (~isfield(shape,'VertexSize'))
        shape.VertexSize = 1;
    end
    if (~isfield(shape,'VertexColor'))
        shape.VertexColor = shape.Color;
    end
    if (~isfield(shape,'LineFlag'))
        shape.LineFlag = 1;
    end
    if (~isfield(shape,'LineColor'))
        shape.LineColor = shape.Color;
    end
end

if strcmp(format,'grfx')
    AddGRFXShape(fid,shape);
elseif strcmp(format,'rvl')
    AddRVLShape(fid,shape);
elseif strcmp(format,'kml')
    AddKMLShape(fid,shape);
end


function AddGRFXShape(fid,shape)
    if (strcmp(shape.ShapeType,'POINT') == 1)
        %write point shape to GRFX file
        fprintf(fid,'<GrGmPoint grType="2000" m_point_symbol_type="0" m_size="%d">\n',shape.VertexSize);
        fprintf(fid,'<GrGmGraphic grType="2000" m_always_save_in_chips="false">\n');
        fprintf(fid,'<GrGraphic grType="2000" m_graphic_name="%s" m_movable="true" m_selectable="true" m_single_point_anchored="true" m_vertices_visible="true" m_visible="true">\n',...
                shape.Name);
        fprintf(fid,'<MthPoint3D name="m_single_point_anchor" x="%f" y="%f" z="0.0"/>\n',shape.Lon(1),shape.Lat(1));
        fprintf(fid,'<PrmPrimitiveList name="m_primitives">\n');
        fprintf(fid,'<PrmEllipse prmType="0" m_angle="0">\n');
        fprintf(fid,'<PrmPrimitive m_scene_space="false" m_is_decoration="false" m_texture_coordinates_exist="false" m_visible="true">\n');
        fprintf(fid,'<MthPoint3DVector name="m_texture_coordinates"/>\n');
        fprintf(fid,'<PrmFaceProperty m_fill_pattern="1" m_face_translucence="1" name="m_face_property">\n');
        fprintf(fid,'<PrmProperty/>\n');
        RGBColor = GetRGBColor(shape.Color);
        fprintf(fid,'<StdRgbColor name="m_face_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',RGBColor(1),...
                     RGBColor(2),RGBColor(3));
        fprintf(fid,'</PrmFaceProperty>\n');
        fprintf(fid,'</PrmPrimitive>\n');
        fprintf(fid,'<MthPoint3D name="m_center" x="0" y="0" z="0"/>\n');
        fprintf(fid,'<MthPoint3D name="m_x_radius_point" x="4" y="0" z="0"/>\n');
        fprintf(fid,'<MthPoint3D name="m_y_radius_point" x="0" y="4" z="0"/>\n');
        fprintf(fid,'</PrmEllipse>\n');
        fprintf(fid,'</PrmPrimitiveList>\n');
        fprintf(fid,'</GrGraphic>\n');
        fprintf(fid,'</GrGmGraphic>\n');
        RGBColor = GetRGBColor(shape.VertexColor);
        fprintf(fid,'<StdRgbColor name="m_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'</GrGmPoint>\n\n');        
    elseif (strcmp(shape.ShapeType,'LINE') == 1)
        fprintf(fid,'<GrGmPolyShape grType="2006" m_closed="false">\n');
        fprintf(fid,'<GrGmGraphic grType="2006" m_always_save_in_chips="false">\n');
        fprintf(fid,'<GrGraphic grType="2006" m_graphic_name="%s">\n',shape.Name);
        RGBColor = GetRGBColor(shape.Color);
        fprintf(fid,'<StdRgbColor name="m_connecting_line_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
                RGBColor(1),RGBColor(2),RGBColor(3));
        if (shape.VertexFlag == 1)    
            fprintf(fid,'<PrmVertexProperty m_vertex_size="%d" m_vertex_symbol="1" name="m_vertex_property">\n',...
                    shape.VertexSize);
            fprintf(fid,'<PrmProperty/>\n');
            RGBColor = GetRGBColor(shape.VertexColor);
            fprintf(fid,'<StdRgbColor name="m_vertex_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
                RBGColor(1),RGBColor(2),RGBColor(3));
            fprintf(fid,'</PrmVertexProperty>\n');
        end
        fprintf(fid,'<PrmEdgeProperty m_edge_thickness="%d" m_edge_pattern="0" m_edge_pattern_changed="true" m_start_vertex_style="0" m_start_vertex_size="15" m_end_vertex_style="0" m_end_vertex_size="15" name="m_edge_property">\n',...
            shape.LineThickness);
        fprintf(fid,'<PrmProperty/>\n');
        RGBColor = GetRGBColor(shape.Color);
        fprintf(fid,'<StdRgbColor name="m_edge_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'<StdRgbColor name="m_start_vertex_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'<StdRgbColor name="m_end_vertex_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'</PrmEdgeProperty>\n');
        fprintf(fid,'<PrmPrimitiveList name="m_primitives">\n');
        fprintf(fid,'<PrmPolyline prmType="4">\n');
        fprintf(fid,'<PrmPolyShape prmType="4" m_vertices_editable="true">\n');
        fprintf(fid,'<PrmPrimitive m_scene_space="false" m_is_decoration="false" m_texture_coordinates_exist="false" m_visible="true">\n');
        fprintf(fid,'<MthPoint3DVector name="m_texture_coordinates"/>\n');
        fprintf(fid,'</PrmPrimitive>\n');
        fprintf(fid,'<MthPoint3DVector name="m_points">\n');
        for i=1:numel(shape.Lat)
            fprintf(fid,'<MthPoint3D name="p%d" x="%f" y="%f" z="0"/>\n',i-1,...
                shape.Lon(i),shape.Lat(i));
        end
        fprintf(fid,'</MthPoint3DVector>\n');
        fprintf(fid,'<StdRgbColorVector name="m_vertex_colors"/>\n');
        fprintf(fid,'</PrmPolyShape>\n');
        fprintf(fid,'</PrmPolyline>\n');
        fprintf(fid,'</PrmPrimitiveList>\n');
        fprintf(fid,'</GrGraphic>\n');
        fprintf(fid,'</GrGmGraphic>\n');
        fprintf(fid,'</GrGmPolyShape>\n');
    elseif (strcmp(shape.ShapeType,'POLYGON') == 1)
        fprintf(fid,'<GrGmPolyShape grType="2005" m_closed="true">\n');
        fprintf(fid,'<GrGmGraphic grType="2005" m_always_save_in_chips="false">\n');
        fprintf(fid,'<GrGraphic grType="2005" m_graphic_name="%s">\n',shape.Name);
        fprintf(fid,'<MthPoint3D name="m_single_point_anchor" x="0" y="0" z="0"/>\n');
        fprintf(fid,'<MthPoint3D name="m_single_point_reference" x="0" y="0" z="0"/>\n');
        RGBColor = GetRGBColor(shape.LineColor);
        fprintf(fid,'<StdRgbColor name="m_connecting_line_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        if (shape.VertexFlag == 1)    
            fprintf(fid,'<PrmVertexProperty m_vertex_size="%d" m_vertex_symbol="1" name="m_vertex_property">\n',...
                    shape.VertexSize);
            fprintf(fid,'<PrmProperty/>\n');
            RGBColor = GetRGBColor(shape.VertexColor);
            fprintf(fid,'<StdRgbColor name="m_vertex_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
                RGBColor(1),RGBColor(2),RGBColor(3));
            fprintf(fid,'</PrmVertexProperty>\n');
        end
        if (shape.LineFlag == 1)
            fprintf(fid,'<PrmEdgeProperty m_edge_thickness="%d" m_edge_pattern="0" m_edge_pattern_changed="true" m_start_vertex_style="0" m_start_vertex_size="15" m_end_vertex_style="0" m_end_vertex_size="15" name="m_edge_property">\n',...
                 shape.LineThickness);
            fprintf(fid,'<PrmProperty/>\n');
            RGBColor = GetRGBColor(shape.LineColor);
            fprintf(fid,'<StdRgbColor name="m_edge_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
                RGBColor(1),RGBColor(2),RGBColor(3));
            fprintf(fid,'<StdRgbColor name="m_start_vertex_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
                RGBColor(1),RGBColor(2),RGBColor(3));
            fprintf(fid,'<StdRgbColor name="m_end_vertex_color" r="%d" g="%d" b="%d" a="1" s="0"/>\n',...
                RGBColor(1),RGBColor(2),RGBColor(3));
            fprintf(fid,'</PrmEdgeProperty>\n'); 
        end
        if (shape.FillFlag == 1)
            fprintf(fid,'<PrmFaceProperty m_fill_pattern="1" m_face_translucence="%f" name="m_face_property">\n',...
                shape.Transparency);
            fprintf(fid,'<PrmProperty/>\n');
            RGBColor = GetRGBColor(shape.FillColor);
            fprintf(fid,'<StdRgbColor name="m_face_color" r="%d" g="%d" b="%d" a="1" s="0"/>',...
                 RGBColor(1),RGBColor(2),RGBColor(3));
            fprintf(fid,'</PrmFaceProperty>\n');
        end
        fprintf(fid,'<PrmPrimitiveList name="m_primitives">\n');
        fprintf(fid,'<PrmPolygon prmType=\"5\">\n');
        fprintf(fid,'<PrmPolyShape prmType="5" m_vertices_editable="true">\n');
        fprintf(fid,'<PrmPrimitive m_scene_space="false" m_is_decoration="false" m_texture_coordinates_exist="false" m_visible="true">\n');
        fprintf(fid,'<MthPoint3DVector name="m_texture_coordinates"/>\n');
        fprintf(fid,'</PrmPrimitive>\n');
        fprintf(fid,'<MthPoint3DVector name="m_points">\n');
        for i=1:numel(shape.Lat)
            fprintf(fid,'<MthPoint3D name="p%d" x="%f" y="%f" z="0"/>\n',i-1,...
                shape.Lon(i),shape.Lat(i));
        end
        fprintf(fid,'</MthPoint3DVector>\n');
        fprintf(fid,'<StdRgbColorVector name=\"m_vertex_colors\"/>\n');
        fprintf(fid,'</PrmPolyShape>\n');
        fprintf(fid,'</PrmPolygon>\n');
        fprintf(fid,'</PrmPrimitiveList>\n');
        fprintf(fid,'</GrGraphic>\n');
        fprintf(fid,'</GrGmGraphic>\n');
        fprintf(fid,'</GrGmPolyShape>\n\n');
    else
        disp 'AddGRFXShape: Unknown Shape Type';
    end
end
    
function AddRVLShape(fid,shape)
    if (strcmp(shape.ShapeType,'POINT') == 1)
        fprintf(fid,'<<\n');
        fprintf(fid,'ClassName	ClassMarker\n');
        fprintf(fid,'Stroke?	true\n');
        RGBColor = GetRGBColor(shape.Color);
        fprintf(fid,'StrokeColor	<-color{%d,%d,%d,0.469047}->\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'ExportShape	<<\n');
        fprintf(fid,'Matrix	[ \n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,']\n');
        fprintf(fid,'Shape	<-GeoCoord{<<\n');
        fprintf(fid,'HorzCoord	<<\n');
        fprintf(fid,'HDatum	WGS-84\n');
        LatStr = GetDMS(shape.Lat(1),'lat');
        LonStr = GetDMS(shape.Lon(1),'lon');
        fprintf(fid,'Lat	<-latitude{%s}-\\>\n',LatStr);
        fprintf(fid,'Lon	<-longitude{%s}-\\>\n',LonStr);
        fprintf(fid,'HUnits	meters\n');
        fprintf(fid,'>>\n');
        fprintf(fid,'>>}->\n');
        fprintf(fid,'Space	LatLon\n');
        fprintf(fid,'InternalSpaceIsLatLon	0\n');
        fprintf(fid,'InternalLineWidth	1 \n');
        fprintf(fid,'InternalLineWidthScale	1.0\n');
        fprintf(fid,'>>\n');
        fprintf(fid,'>>\n\n');
    elseif (strcmp(shape.ShapeType,'LINE') == 1)
        fprintf(fid,'<<\n');
        fprintf(fid,'ClassName	ClassPolyLineTool\n');
        fprintf(fid,'ExplorerTitle	PolyLineTool\n');
        fprintf(fid,'VisibilityType	1\n');
        fprintf(fid,'Stroke?	true\n');
        RGBColor = GetRGBColor(shape.Color);
        fprintf(fid,'StrokeColor	<-color{%d,%d,%d,0.469047}->\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'Fill?	false\n');
        fprintf(fid,'ExportShape	<<\n');
        fprintf(fid,'Matrix	[ \n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,']\n');
        fprintf(fid,'Shape	[\n');
        for i=1:numel(shape.Lat)
            fprintf(fid,'<-GeoCoord{<<\n');
            fprintf(fid,'HorzCoord	<<\n');
            fprintf(fid,'HDatum	WGS-84\n');
            LatStr = GetDMS(shape.Lat(i),'lat');
            LonStr = GetDMS(shape.Lon(i),'lon');
            fprintf(fid,'Lat	<-latitude{%s}-\\>\n',LatStr);
            fprintf(fid,'Lon	<-longitude{%s}-\\>\n',LonStr);
            fprintf(fid,'HUnits	meters\n');
            fprintf(fid,'>>\n');
            fprintf(fid,'>>}->\n');
        end
        fprintf(fid,']\n');
        fprintf(fid,'Space	LatLon\n');
        fprintf(fid,'InternalMatrix	[ \n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,']\n');
        fprintf(fid,'InternalSpaceIsLatLon	0\n');
        fprintf(fid,'InternalLineWidth	%d\n',shape.LineThickness);
        fprintf(fid,'InternalLineWidthScale	1.0\n');
        fprintf(fid,'>>\n');
        fprintf(fid,'>>\n\n');
    elseif (strcmp(shape.ShapeType,'POLYGON') == 1)
        fprintf(fid,'<<\n');
        fprintf(fid,'ClassName	ClassClosedPolygonTool\n');
        fprintf(fid,'VisibilityType	1\n');
        fprintf(fid,'Stroke?	true\n');
        RGBColor = GetRGBColor(shape.Color);
        fprintf(fid,'StrokeColor	<-color{%d,%d,%d,0.469047}->\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'Fill?	%d\n',shape.FillFlag);
        RGBColor = GetRGBColor(shape.FillColor);
        fprintf(fid,'FillColor	<-color{%d,%d,%d,0.469047}->\n',...
            RGBColor(1),RGBColor(2),RGBColor(3));
        fprintf(fid,'FillPattern	1\n');
        fprintf(fid,'FillOpacity	false\n');
        fprintf(fid,'LineJoin	1\n');
        fprintf(fid,'ExportShape	<<\n');
        fprintf(fid,'Matrix	[ \n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'1.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,'0.0\n');
        fprintf(fid,']\n');
        fprintf(fid,'Shape	[');
        for i=1:numel(shape.Lat)
            fprintf(fid,'<-GeoCoord{<<\n');
            fprintf(fid,'HorzCoord	<<\n');
            fprintf(fid,'HDatum	WGS-84\n');
            LatStr = GetDMS(shape.Lat(i),'lat');
            LonStr = GetDMS(shape.Lon(i),'lon');
            fprintf(fid,'Lat	<-latitude{%s}-\\>\n',LatStr);
            fprintf(fid,'Lon	<-longitude{%s}-\\>\n',LonStr);
            fprintf(fid,'HUnits	meters\n');
            fprintf(fid,'>>\n');
            fprintf(fid,'>>}->\n');
        end
        fprintf(fid,']\n');
        fprintf(fid,'Space	LatLon\n');
        fprintf(fid,'InternalSpaceIsLatLon	0\n');
        fprintf(fid,'InternalLineWidth	%d\n',shape.LineThickness);
        fprintf(fid,'InternalLineWidthScale	1.0\n');
        fprintf(fid,'>>\n');
        fprintf(fid,'>>\n\n');
    else
        disp 'AddRVLShape: Unknown Shape Type';
    end
end

function AddKMLShape(fid,shape)
    if (strcmp(shape.ShapeType,'POINT') == 1)
        fprintf(fid,'<Placemark>\n');
        fprintf(fid,'<name>%s</name>\n',shape.Name);
        fprintf(fid,'<Style>\n');
        fprintf(fid,'<IconStyle>\n');
        fprintf(fid,'<Icon>\n');
        fprintf(fid,'  <href>%s</href>\n', GetKMLPushpinFile(shape.Color));
        fprintf(fid,'  <scale>1.0</scale>\n');
        fprintf(fid,'</Icon>\n');
        fprintf(fid,'</IconStyle>\n');
        fprintf(fid,'</Style>\n');
        if isfield(shape,'Date') %test to see if user set the Date field
            fprintf(fid,'<TimeStamp><when>%s</when></TimeStamp>\n',shape.Date);
        end
        if isfield(shape,'StartTime') && isfield(shape,'StopTime')
            %specify time span
            fprintf(fid,'<TimeSpan>\n<begin>%s</begin>\n',shape.StartTime);
            fprintf(fid,'<end>%s</end>\n</TimeSpan>\n',shape.StopTime);
        end
        fprintf(fid,'<description>%s</description>\n',shape.Description);
        fprintf(fid,'<Point>\n');
        if AltFlag
            fprintf(fid,'<altitudeMode>absolute</altitudeMode>\n');
        else
            fprintf(fid,'<altitudeMode>clampToGround</altitudeMode>\n');
        end
        fprintf(fid,'<coordinates>%f,%f,%f</coordinates>\n',...
            shape.Lon(1),shape.Lat(1),shape.Alt(1));
        fprintf(fid,'</Point>\n');
        fprintf(fid,'</Placemark>\n\n');        
    elseif (strcmp(shape.ShapeType,'LINE') == 1) 
        fprintf(fid,'<Placemark>\n');
        fprintf(fid,'<name>%s</name>\n',shape.Name);
        fprintf(fid,'<Style>\n');
        fprintf(fid,'<LineStyle>\n');
        fprintf(fid,'<width>%d</width>\n',shape.LineThickness);
        Color = GetKMLColor(shape.Color,1);
        fprintf(fid,'<color>%s</color>\n',Color);
        fprintf(fid,'</LineStyle>\n');
        fprintf(fid,'</Style>\n');        
        if isfield(shape,'Date') %test to see if user set the Date field
            fprintf(fid,'<TimeStamp><when>%s</when></TimeStamp>\n',shape.Date);
        end
        if isfield(shape,'StartTime') && isfield(shape,'StopTime')
            %specify time span
            fprintf(fid,'<TimeSpan>\n<begin>%s</begin>\n',shape.StartTime);
            fprintf(fid,'<end>%s</end>\n</TimeSpan>\n',shape.StopTime);
        end
        fprintf(fid,'<description>%s</description>\n',shape.Description);
        fprintf(fid,'<LineString>\n');
        if AltFlag
            fprintf(fid,'<altitudeMode>absolute</altitudeMode>\n');
        else
            fprintf(fid,'<altitudeMode>clampToGround</altitudeMode>\n');
        end
        fprintf(fid,'<coordinates>\n');
        for i=1:numel(shape.Lat)
            fprintf(fid,'%f,%f,%f\n',shape.Lon(i),shape.Lat(i),shape.Alt(i));
        end
        fprintf(fid,'</coordinates>\n');
        fprintf(fid,'</LineString>\n');
        fprintf(fid,'</Placemark>\n\n');        
    elseif (strcmp(shape.ShapeType,'POLYGON') == 1)
        fprintf(fid,'<Placemark>\n');
        fprintf(fid,'<name>%s</name>\n',shape.Name);
        fprintf(fid,'<Style>\n');
        fprintf(fid,'<LineStyle>\n');
        fprintf(fid,'<width>%d</width>\n',shape.LineThickness);
        Color = GetKMLColor(shape.Color,1);
        fprintf(fid,'<color>%s</color>\n',Color);
        fprintf(fid,'</LineStyle>\n');
        fprintf(fid,'<PolyStyle>\n');
        fprintf(fid,'<width>%d</width>\n',shape.LineThickness);
        Color = GetKMLColor(shape.FillColor,shape.Transparency);
        fprintf(fid,'<color>%s</color>\n',Color);
        fprintf(fid,'</PolyStyle>\n');
        fprintf(fid,'</Style>\n');
        if isfield(shape,'Date') %test to see if user set the Date field
            fprintf(fid,'<TimeStamp><when>%s</when></TimeStamp>\n',shape.Date);
        end
        if isfield(shape,'StartTime') && isfield(shape,'StopTime')
            %specify time span
            fprintf(fid,'<TimeSpan>\n<begin>%s</begin>\n',shape.StartTime);
            fprintf(fid,'<end>%s</end>\n</TimeSpan>\n',shape.StopTime);
        end
        fprintf(fid,'<description>%s</description>\n',shape.Description);
        fprintf(fid,'<Polygon>\n');
        %%% These cause bad ground clipping effects in Google Earth
        %%% 4.0.2737
        %%% especially when zooming.  Needed?
        %%%fprintf(fid,'<extrude>1</extrude>\n');
        if AltFlag
            fprintf(fid,'<altitudeMode>absolute</altitudeMode>\n');
        else
            fprintf(fid,'<altitudeMode>clampToGround</altitudeMode>\n');
        end
        fprintf(fid,'<outerBoundaryIs>\n');
        fprintf(fid,'<LinearRing>\n');
        fprintf(fid,'<coordinates>\n');
        for i=1:numel(shape.Lat)
            fprintf(fid,'%f,%f,%f\n',shape.Lon(i),shape.Lat(i),shape.Alt(i));
        end
        fprintf(fid,'</coordinates>\n');
        fprintf(fid,'</LinearRing>\n');
        fprintf(fid,'</outerBoundaryIs>\n');
        fprintf(fid,'</Polygon>\n');
        fprintf(fid,'</Placemark>\n\n');
    elseif (strcmp(shape.ShapeType,'OVERLAY') == 1)
        fprintf(fid,'<GroundOverlay>\n');
        fprintf(fid,'<name>%s</name>\n',shape.Name);
        if isfield(shape,'Date') %test to see if user set the Date field
            fprintf(fid,'<TimeStamp><when>%s</when></TimeStamp>\n',shape.Date);
        end
        if isfield(shape,'StartTime') && isfield(shape,'StopTime')
            %specify time span
            fprintf(fid,'<TimeSpan>\n<begin>%s</begin>\n',shape.StartTime);
            fprintf(fid,'<end>%s</end>\n</TimeSpan>\n',shape.StopTime);
        end
        fprintf(fid,'<Icon>\n');
        fprintf(fid,'<href>%s</href>\n',shape.OverlayFile);
        fprintf(fid,'</Icon>\n');
        fprintf(fid,'<LatLonBox>\n');
        fprintf(fid,'<north>%f</north>\n',shape.OverlayBox(1));
        fprintf(fid,'<south>%f</south>\n',shape.OverlayBox(2));
        fprintf(fid,'<east>%f</east>\n',shape.OverlayBox(3));
        fprintf(fid,'<west>%f</west>\n',shape.OverlayBox(4));
        fprintf(fid,'</LatLonBox>\n');
        fprintf(fid,'</GroundOverlay>\n\n');
    else
        disp 'AddKMLShape: Unknown Shape Type';
    end
end

end


function RGBColor = GetRGBColor(Color)
    %if RGB Triplet is passed in, then just return it
    if (Color(1) >= 0 && Color(1) <= 1 &&...
         Color(2) >= 0 && Color(2) <= 1 &&...
         Color(3) >= 0 && Color(3) <= 1)
        RGBColor(1) = Color(1);
        RGBColor(2) = Color(2);
        RGBColor(3) = Color(3);  
        return;
    end
    
    switch lower(Color)
        case 'red'
            RGBColor(1) = 1;
            RGBColor(2) = 0;
            RGBColor(3) = 0;
        case 'orange'
            RGBColor(1) = 1;
            RGBColor(2) = 0.4;
            RGBColor(3) = 0;
        case 'yellow'
            RGBColor(1) = 1;
            RGBColor(2) = 1;
            RGBColor(3) = 0;
        case 'green'
            RGBColor(1) = 0;
            RGBColor(2) = 1;
            RGBColor(3) = 0;
        case 'cyan'
            RGBColor(1) = 0;
            RGBColor(2) = 1;
            RGBColor(3) = 1;
        case 'blue'
            RGBColor(1) = 0;
            RGBColor(2) = 0;
            RGBColor(3) = 1;
        case 'magenta'
            RGBColor(1) = 1;
            RGBColor(2) = 0;
            RGBColor(3) = 1;
        case 'white'
            RGBColor(1) = 1;
            RGBColor(2) = 1;
            RGBColor(3) = 1;
        case 'black'
            RGBColor(1) = 1;
            RGBColor(2) = 1;
            RGBColor(3) = 1;    
        otherwise                      
            %set red as default
            RGBColor(1) = 1;
            RGBColor(2) = 0;
            RGBColor(3) = 0;      
    end
return; 
end

function out_string = GetDMS(Deg,LatLon)
    out_string = latlonstr(Deg, LatLon, 'num_units', 3, ...
        'delimiter', {':',':',' '}, 'include_symbols', false, ...
        'signed', false, 'precision', 2);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////