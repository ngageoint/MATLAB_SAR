function process_phd(hObject,handles)
% Condition PHD as specified (center, deweight, left/right flip)
for ii = 1:size(handles.complex_data,3) %treat PHD for each polarization separately
    cdata = handles.complex_data(:,:,ii);
    meta = handles.meta;
    if all(handles.manual_offset==0)
        % If possible deweight in an already deskewed dimension before deskewing in
        % the other dimension.  We will deweight the other dimension after deskew.
        if all(meta.Grid.Col.DeltaKCOAPoly(:) == 0) && ...
                all(meta.Grid.Row.DeltaKCOAPoly(:) == 0) && ...
                get(handles.UniformWeighting,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
            cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
        elseif all(meta.Grid.Col.DeltaKCOAPoly(:) == 0) && ...
                get(handles.UniformWeighting,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
        elseif all(meta.Grid.Row.DeltaKCOAPoly(:) == 0) && ...
                get(handles.UniformWeighting,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
        end
        % Frequency support deskew
        if get(handles.DeskewSlow,'Value')
            [ DeltaKCOAPoly, az_coords_m, rg_coords_m, fft_sgn ] = deskewparams( meta, 1 );
            if any(DeltaKCOAPoly(:)~=0) % Do we need to deskew frequency support?
                [cdata, new_DeltaKCOAPoly] = deskewmem(cdata, DeltaKCOAPoly, az_coords_m, rg_coords_m, 1, fft_sgn);
                % Deskewing might shift frequency support in other direction.
                if get(handles.UniformWeighting,'Value')
                    cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
                end
            elseif ~get(handles.DeskewFast,'Value') % If no deskew, at least attempt to recenter
                new_DeltaKCOAPoly = meta.Grid.Row.DeltaKCOAPoly;
            else
                new_DeltaKCOAPoly = 0;
            end
            if any(new_DeltaKCOAPoly(:)~=0)
                % In general, we cannot deskew in both directions at once.
                % However, we can at least recenter the non-deskewed dimension
                % with a simple shift (not a proper deskew, which is a spatially
                % variant shift).
                deltaKCOA = sicd_polyval2d(new_DeltaKCOAPoly,...
                    az_coords_m(round(numel(az_coords_m)/2)),...
                    rg_coords_m(round(numel(rg_coords_m)/2))); % Get shift at center of AOI
                % deltaKCOA is scalar so it results in a spatially invariant
                % shift which will not affect the deskew diresction.
                cdata = deskewmem(cdata, deltaKCOA, az_coords_m, rg_coords_m, 2, meta.Grid.Row.Sgn);
            end
        end
        if get(handles.DeskewFast,'Value')
            [ DeltaKCOAPoly, az_coords_m, rg_coords_m, fft_sgn ] = deskewparams( meta, 2 );
            if any(DeltaKCOAPoly(:)~=0) % Do we need to center frequency support?
                [cdata, new_DeltaKCOAPoly] = deskewmem(cdata, DeltaKCOAPoly, az_coords_m, rg_coords_m, 2, fft_sgn);
                if get(handles.UniformWeighting,'Value')
                    cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
                end
            elseif ~get(handles.DeskewSlow,'Value') % If no deskew, at least attempt to recenter
                new_DeltaKCOAPoly = meta.Grid.Col.DeltaKCOAPoly;
            else
                new_DeltaKCOAPoly = 0;
            end
            if any(new_DeltaKCOAPoly(:)~=0)
                % See comments above for deskew and the resulting shift
                deltaKCOA = sicd_polyval2d(new_DeltaKCOAPoly,...
                    az_coords_m(round(numel(az_coords_m)/2)),...
                    rg_coords_m(round(numel(rg_coords_m)/2))); % Get shift at center of AOI
                cdata = deskewmem(cdata, deltaKCOA, az_coords_m, rg_coords_m, 1, meta.Grid.Col.Sgn);
            end
        end
    else
        cdata = handles.fft_im(circshift(handles.fft_sp(cdata),handles.manual_offset)); % Manual shift
        if get(handles.UniformWeighting,'Value')
            cdata = deweightmem(cdata, handles.weight_fun_az, handles.AzPad, 1);
            cdata = deweightmem(cdata, handles.weight_fun_rng, handles.RnPad, 2);
        end
    end
    
    phd(:,:,ii) = fftshift(handles.fft_sp(cdata.'));
    try
        [temp,k_a,k_r] = pfa2inv(phd(:,:,ii),meta);
        phd_inv(:,:,ii) = double(temp);
        handles.k_a = k_a;
        handles.k_r = k_r;
        handles.phasehistory_inv = phd_inv;
        set(handles.InversePolar,'Enable','on');
    catch        
        set(handles.InversePolar,'Enable','off');
        set(handles.InversePolar,'Value',0);
    end
end
% Uniform weighting
if get(handles.UniformWeighting,'Value') && ...
        isfield(meta.Grid.Row,'ImpRespBW') && isfield(meta.Grid.Col,'ImpRespBW')
    % Adjust impulse response width metadata for new weighting
    meta.Grid.Col.ImpRespWid = .886/meta.Grid.Col.ImpRespBW;
    meta.Grid.Row.ImpRespWid = .886/meta.Grid.Row.ImpRespBW;
end

handles.phasehistory = phd;

%store chip.  We will need to fix the contrast in resolution mode, so when
%it is selected, we will store the displayed sigmas above mean
% for ii = 1:size(phd,3)
%     chipmag(:,:,ii) = abs(handles.fft_im(phd(:,:,ii)));
% end
[chipmag,handles] = ApToolmakeDisplayable(handles,phd);
handles.chip = chipmag;

%Draw Image 
remaps = get(handles.RemapCombo,'String');
remap = remaps{get(handles.RemapCombo,'Value')};
if strcmp(remap,'linearremap')
    chipmean = mean(chipmag(:));
    chipstd = std(chipmag(:));
    handles.apiSP1.replaceImage(chipmag,...
        'DisplayRange', [0 chipmean+3*chipstd], 'PreserveView', true);
else
    handles.apiSP1.replaceImage(chipmag, 'PreserveView', true);
end

%Draw PHD
if isfield(handles,'polar_PHD_idx')
    if get(handles.InversePolar,'Value')
        phdmag = abs(phd_inv(:,:,handles.polar_PHD_idx));
    else
        phdmag = abs(phd(:,:,handles.polar_PHD_idx));
    end
else
    if get(handles.InversePolar,'Value')
        phdmag = abs(phd_inv);
    else
        phdmag = abs(phd);
    end
end
phdmag(isnan(phdmag)) = 0;
phdmean = mean(phdmag(:));
phdmag = 10.*phdmag./phdmean;
phdmean = 10;
phdstd = std(phdmag(:));

%if flight is left then the phd needs to be flipped horizontally to
%represent time going right (which we do in this country...)
if isfield(meta,'SCPCOA') && isfield(meta.SCPCOA,'SideOfTrack') && ...
        strcmp(meta.SCPCOA.SideOfTrack,'L')
    phdmag = fliplr(phdmag);
    flipXaxis = 0;
else
    flipXaxis = 1;
end

%make axes labels
try
    if isfield(handles,'phdYaxis')
        handles = rmfield(handles,'phdYaxis');
    end
    freq_width = (1/handles.meta.Grid.Row.SS)*(SPEED_OF_LIGHT/2);
    freq_ctr = handles.meta.Grid.Row.KCtr*(SPEED_OF_LIGHT/2);
    freq_limits = freq_ctr + ([-1 1]*freq_width/2);
    if isfield(handles.meta,'PFA') && isfield(handles.meta.PFA,'SpatialFreqSFPoly')
        freq_limits = freq_limits/handles.meta.PFA.SpatialFreqSFPoly(1);
    end
    if freq_ctr > 1e9
        freq_limits = freq_limits/1e9;
        handles.freqUnits = 'GHz';
    else
        freq_limits = freq_limits/1e6;
        handles.freqUnits = 'MHz';
    end
    handles.phdYaxis = linspace(freq_limits(2), freq_limits(1), size(phdmag,1));
end
try
    if isfield(handles,'phdXaxis')
        handles = rmfield(handles,'phdXaxis');
    end
    angle_width = (1/handles.meta.Grid.Col.SS) / handles.meta.Grid.Row.KCtr;
    if isfield(handles.meta.Grid.Col,'KCtr')
        angle_ctr = handles.meta.Grid.Col.KCtr;
    else
        angle_ctr = 0;
    end
    angle_limits = angle_ctr+([-1 1]*angle_width/2);
    if flipXaxis
        angle_limits = angle_limits([2 1]);
    end
    handles.phdXaxis = atand(linspace(angle_limits(1), angle_limits(2), size(phdmag,2)));
end

set(handles.phd,'XLim',[1 size(phdmag,2)]);
set(handles.phd,'YLim',[1 size(phdmag,1)]);

%update PHD image
if ~isfield(handles,'phdhandle')
   handles.phdhandle = imagesc(phdmag,'parent',handles.phd);   
else
    set(handles.phdhandle,'cdata',phdmag)
end
set(handles.phd,'Colormap',jet)

handles.phdmagSize = size(phdmag);

ApToolupdatePHD_axesLabels(handles)

% Update handles structure
guidata(hObject, handles);
