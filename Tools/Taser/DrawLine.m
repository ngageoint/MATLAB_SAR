function [ Distance Bearing Line ] = DrawLine( obj, Line )

if Line ~= []
    delete(Line);
    Line = [];
end
Line = imline(obj.mitm_hand.AxesHandle);
setColor(Line,'Cyan');
kids = get(Line,'Children');
set(kids,'LineWidth',2);
set(kids,'MarkerSize',1);
%set(obj.DrawLine,'String','Remove Line');
%add callback for user moving line
api_handle = iptgetapi(Line);
api_handle.addNewPositionCallback(@newlinepos);
pos = getPosition(Line);
Distance = ComputeDistance(obj.mitm_hand, pos);
meta = obj.mitm_hand.Metadata{obj.mitm_hand.Frame};
Bearing = 90 + atan2((pos(2,2)-pos(1,2))*meta.Grid.Row.SS,(pos(2,1)-pos(1,1))*meta.Grid.Row.SS)*180/pi;
if Bearing < 0
    Bearing = 360 + Bearing;
end
if Bearing >= 360
    Bearing = Bearing -360;
end

function newlinepos(pos)

handles = guidata(gcbo);
distance = ComputeDistance(handles.mitm_hand, pos);
meta = handles.mitm_hand.Metadata{handles.mitm_hand.Frame};
Bearing = 90 + atan2((pos(2,2)-pos(1,2))*meta.Grid.Row.SS,(pos(2,1)-pos(1,1))*meta.Grid.Row.SS)*180/pi;
if Bearing < 0
    Bearing = 360 + Bearing;
end
if Bearing >= 360
    Bearing = Bearing -360;
end
set(handles.txtMeasure, 'String', [num2str(distance) ' m, ' num2str(Bearing) ' degrees']);