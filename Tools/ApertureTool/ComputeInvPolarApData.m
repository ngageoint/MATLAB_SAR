function Ap = ComputeInvPolarApData(XBounds,YBounds,handles)
%COMPUTEINVPOLARAPDATA Filter PHD from k_r/k_a defined AOI

[ny, nx, nz] = size(handles.phasehistory);
Ap = zeros(ny,nx,nz);

%convention is upside down...
YBounds = [ny-YBounds(2) ny-YBounds(1)];

%define k_r/k_a for selected AOI coordinates
k_a = linspace(handles.k_a(1),handles.k_a(2),nx);
k_r = linspace(handles.k_r(1),handles.k_r(2),ny);

%construct rectangle with 10 points per side in k_a, k_r space
ka(1:10) = k_a(XBounds(1));
kr(1:10) = linspace(k_r(YBounds(1)),k_r(YBounds(2)),10);
ka(11:20) = linspace(k_a(XBounds(1)),k_a(XBounds(2)),10);
kr(11:20) = k_r(YBounds(2));
ka(21:30) = k_a(XBounds(2));
kr(21:30) = linspace(k_r(YBounds(2)),k_r(YBounds(1)),10);
ka(31:40) = linspace(k_a(XBounds(2)),k_a(XBounds(1)),10);
kr(31:40) = k_r(YBounds(1));

%compute k_u and k_v spacing from metadata
center_freq = handles.meta.Grid.Row.KCtr;
ss = [handles.meta.Grid.Row.SS handles.meta.Grid.Col.SS];
k_v = center_freq + ((1./ss(2)) * [-1 1] / 2);
k_u = (1./ss(1)) * [-1 1] / 2;
k_v = linspace(k_v(1), k_v(2), ny)';
k_u = linspace(k_u(1), k_u(2), nx)';

%loop through points and convert to pfa space
xpos = zeros(length(ka),1);
ypos = zeros(length(kr),1);
for i=1:length(ka)
    ku_temp = sind(ka(i))*kr(i);
    kv_temp = cosd(ka(i))*kr(i);
    %get pixel position in ku/kv space
    [~,index1] = min(abs(k_u-ku_temp));
    [~,index2] = min(abs(k_v-kv_temp));
    xpos(i) = index1;
    ypos(i) = ny - index2;    
end

%now convert polygon to mask
BW = poly2mask(xpos,ypos,ny,nx);

%invert mask if inverse selection
if get(handles.InverseSelection,'Value')
    BW = -1*(BW-1);
end

%filter PHD by BW
for i=1:nz
    Ap(:,:,i) = handles.phasehistory(:,:,i).*BW;
end

end

