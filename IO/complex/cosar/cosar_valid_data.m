function [ valid_data ] = cosar_valid_data( filename, az_symmetry )
%COSAR_VALID_DATA Extract valid data region from TerraSAR-X COSAR file format
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('az_symmetry','var')
    az_symmetry = false; % Most TSX data is right-looking
end

% Open file
fid=fopen(filename,'r','b');
fseek(fid,8,'bof');
rs=fread(fid,1,'int32'); % Range samples
az=fread(fid,1,'int32'); % Azimuth samples

% Look for valid data in azimuth extents
fseek(fid,(rs+2)*4,'bof'); % Skip burst annotation
az_annotation=fread(fid,[(rs+2) 3],'int32'); % Valid data extents
az_min=az_annotation(3:end,2).'; % Trim off filler values
az_max=az_annotation(3:end,3).'; % Trim off filler values

% Look for valid data is range extents
fseek(fid,((rs+2)*4*4),'bof'); % Skip burst annotation
rg_annotation=fread(fid,[2 az],'2*int32', rs*4); % Valid data extents
rg_min = rg_annotation(1,:);
rg_max = rg_annotation(2,:);

% Close file
fclose(fid);

% COSAR defines data validity with bounds for each row and column.  For
% data to be valid it must be valid in both row and column.  Here we throw
% out points that are outside the bounding polygon.
rows = 1:rs;
azmin_valid = rows>=rg_min(az_min(rows)) & rows<=rg_max(az_min(rows));
azmax_valid = rows>=rg_min(az_max(rows)) & rows<=rg_max(az_max(rows));
cols = 1:az;
rgmin_valid = cols>=az_min(rg_min(cols)) & cols<=az_max(rg_min(cols));
rgmax_valid = cols>=az_min(rg_max(cols)) & cols<=az_max(rg_max(cols));
% Define a countour of valid points that surrounds the data.
boundary_col = [az_min(azmin_valid) cols(rgmin_valid) az_max(azmax_valid) cols(rgmax_valid)];
if az_symmetry
    boundary_col = az - boundary_col + 1;
end
boundary_row = [rows(azmin_valid) rg_min(rgmin_valid) rows(azmax_valid) rg_max(rgmax_valid)];
% Create convex hull that contains these boundary points
K = convhull(boundary_col,boundary_row,'simplify',true);

% Reorder convhull output to SICD standard.  convhull already orients
% clockwise, just need to start on correct vertex.
% From SICD spec: Vertex 1 determined by: (1) minimum row index, (2)
% minimum column index if 2 vertices with minimum row index.
K = K(1:(end-1)); % Last output from convhull is repeat of first
minrow = min(boundary_row(K));
indr = find(boundary_row(K)==minrow); % Might be more than one
[mincol, indc] = min(boundary_col(K(indr)));
first_vertex = indr(indc);
K = circshift(K, 1 - first_vertex);

% Return SICD structure
for i = 1:numel(K)
    valid_data.Vertex(i).Row = boundary_row(K(i)) - 1; % SICD indices are zero-based, not one-based like MATLAB
    valid_data.Vertex(i).Col = boundary_col(K(i)) - 1;
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////