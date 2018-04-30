function snail_trail( filename, callback_function )
%SNAIL_TRAIL Allows user to select AOIs repeatedly and remember/display them
%
% Example:
% snail_trail(filename, @disp)
%
% This example just prints bounds of selected areas, but you could do
% processing in the callback as well.
%
% TODO: Make this work for multi-image files.
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Setup layout
fig_hand = figure('MenuBar','none','Toolbar','none');
uip_hand = uipanel('Parent',fig_hand,'Position',[0 0 1 1],'BorderType','none');
mitm_hand = hg_mitm_viewer(uip_hand);
aoiButton = javaObjectEDT('javax.swing.JButton','Select AOI');
set(handle(aoiButton,'callbackproperties'), 'ActionPerformedCallback', ...
    @(obj,eventdata) process_aoi());

% Open file
mitm_hand.openFile(filename, true);
% Might be cleaner if we fit figure shape to image aspect ratio first like in TASER.
mitm_hand.Zoom = 'fit';
mitm_hand.PreChangeViewFcn = @() saveShapesNativeCoords();
mitm_hand.PostChangeViewFcn = @() restoreShapesLocalCoords();

% Add buttons to toolbar
main_toolbar = mitm_hand.main_toolbar;
main_toolbar.addSeparator;
main_toolbar.add(aoiButton);
main_toolbar.repaint;
main_toolbar.revalidate;

% Setup draggable rectangle for AOI selection
datasize = double([mitm_hand.Metadata{mitm_hand.Frame}.ImageData.NumCols ...
            mitm_hand.Metadata{mitm_hand.Frame}.ImageData.NumRows]);
native_coords = [datasize(1)/4 datasize(2)/4 datasize(1)/2-1 datasize(2)/2-1];
init_coords = mitm_hand.native2axescoords(native_coords(1:2));
init_coords(3:4) = native_coords(3:4)/max(mitm_hand.Zoom,1);
aoi_imrect = imrect(mitm_hand.AxesHandle, init_coords);
init_constraint = mitm_hand.native2axescoords(datasize+1);
% Force draggable rectangle to be within image, rather than entire figure/axes
setPositionConstraintFcn(aoi_imrect, makeConstrainToRectFcn('imrect', ...
    [1 init_constraint(1)], [1 init_constraint(2)]));
processed_areas = {};


    % This function shades all processed areas as reddish.
    function data_out = shade_processed_areas(data_in)
        % Doesn't currently handle muli-frame datasets, in which the masks
        % of processed area are almost certainly different for each frame.
        data_out = repmat(data_in, [1 1 3]);
        mask = zeros(size(data_in));
        for i = 1:numel(processed_areas)
            ul = mitm_hand.native2axescoords(processed_areas{i}(1:2));
            lr = ul + processed_areas{i}(3:4)/max(mitm_hand.Zoom,1) - 1;
            mask = mask | poly2mask([ul(1), lr(1), lr(1), ul(1)], ...
                [ul(2), ul(2), lr(2), lr(2)], size(data_in,1), size(data_in,2));
        end
        % Shade processed areas red
        data_out(:,:,2) = data_out(:,:,2) .* ~mask;
        data_out(:,:,3) = data_out(:,:,3) .* ~mask;
    end

    % This function uses the draggable rectangle to start AOI processing.
    function process_aoi
        % Get AOI
        aoi_screen = aoi_imrect.getPosition();
        % Convert to image coordinate space
        ulcorner = round(mitm_hand.axescoords2native(aoi_screen(1:2)));
        lrcorner = round(mitm_hand.axescoords2native(aoi_screen(1:2) + aoi_screen(3:4) - 1));
        % Return values (UL corner and size) to calling function
        callback_args{1}(1:2) = ulcorner;
        callback_args{1}(3:4) = lrcorner-ulcorner+1;
        callback_args{2} = mitm_hand.FrameSubFileIndices{mitm_hand.Frame}; % Which frame was AOI selected in?
        try
            callback_function(callback_args{1});
            processed_areas{end+1} = callback_args{1};  % Do we need to store frame as well?
            % This function never changes (only the processed_areas
            % variable it references), but setting it here, even if to the
            % same function it was before, causes the viewer to actually
            % apply it.
            mitm_hand.DataTransformFcn = @shade_processed_areas;
        end
    end

    % Rember draggable rectangle before zoom/pan
    function saveShapesNativeCoords()
        local_coords = aoi_imrect.getPosition();
        native_coords = mitm_hand.axescoords2native(local_coords(1:2));
        native_coords(3:4) = local_coords(3:4)*max(mitm_hand.Zoom,1);
    end

    % Restore draggable rectangle before zoom/pan
    % This is messier than one would hope because of the MATLAB quirks.
    function restoreShapesLocalCoords()
        oldzoom = get(zoom(fig_hand),'Enable');
        oldpan = get(pan(fig_hand),'Enable');
        delete(aoi_imrect);
        local_coords = mitm_hand.native2axescoords(native_coords(1:2));
        local_coords(3:4) = native_coords(3:4)/max(mitm_hand.Zoom,1);
        aoi_imrect = imrect(mitm_hand.AxesHandle, local_coords);
        ds = double([mitm_hand.Metadata{mitm_hand.Frame}.ImageData.NumCols ...
            mitm_hand.Metadata{mitm_hand.Frame}.ImageData.NumRows]);
        constraint = mitm_hand.native2axescoords(ds+1);
        setPositionConstraintFcn(aoi_imrect, makeConstrainToRectFcn('imrect', ...
            [1 constraint(1)], [1 constraint(2)]));
        % Restore old state which was changed when we redrew the shapes
        try % Doesn't work when initiated from zoom callback.
            zoom(fig_hand,oldzoom);
            pan(fig_hand, oldpan);
        end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////