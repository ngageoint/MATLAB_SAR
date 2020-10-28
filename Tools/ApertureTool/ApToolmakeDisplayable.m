function [out,handles] = ApToolmakeDisplayable(handles,phd)
% Do all processing to convert raw complex data into a real format
% MATLAB can display

for ii = 1:size(phd,3)
    out(:,:,ii) = handles.fft_im(phd(:,:,ii));
end

reps = get(handles.image_rep_combo,'String');
rep = reps{get(handles.image_rep_combo,'Value')};

if ~any(strcmpi(rep,{'amplitude' 'Pauli' 'AlphaEntropy' 'Stokes'}))
    error('ApertureTool:ImageRep','ApertureTool: amplitude, Pauli, AlphaEntropy, and Stokes are the only supported image representations');
end
if size(out,3) > 1
    if strcmpi(rep,'amplitude')
        %use single polarization channel (treat as single band)
        out = out(:,:,handles.polar_PHD_idx);
    else
        switch size(out,3)
            case 1 % Single band image; nothing to do
            case 2 % Dual-pol
                co_index=[find(strcmpi('H:H',handles.TxRcvPolarizationProc))...
                    find(strcmpi('V:V',handles.TxRcvPolarizationProc))];
                if isempty(co_index) %could be RHC or LHC
                    co_index=1;
                    cross_index=2;
                else
                    if length(co_index)>1 % Catch HH/VV
                        cross_index=co_index(2);
                        co_index=co_index(1);
                    else
                        cross_index=[find(strcmpi('H:V',handles.TxRcvPolarizationProc))...
                            find(strcmpi('V:H',handles.TxRcvPolarizationProc))];
                    end
                end
                if strcmpi(rep,'AlphaEntropy')
                    [out(:,:,1),out(:,:,2),out(:,:,3)] =ComputeAlphaEntropy(out(:,:,co_index),...
                        out(:,:,cross_index),[],[]);
                elseif strcmpi(rep,'Pauli')
                    out=cat(3,abs(out(:,:,co_index)),...
                        abs(out(:,:,cross_index)),...
                        abs(out(:,:,co_index)));
                elseif strcmpi(rep,'Stokes')
                    [S,m,psi,chi,del] = ComputeStokesDual(out(:,:,co_index),out(:,:,cross_index),5);
                    Cutoff = 2.5*std(psi(:));
                    psi(psi>Cutoff) = Cutoff;
                    psi(psi<-Cutoff) = -Cutoff;
                    psi = abs(psi);
                    psi = psi./Cutoff;
                    %cut off HSV at 175/255 (Red to Blue)
                    psi = psi.*(175/255);
                    MinM = min(m(:));
                    m = m - MinM;
                    m = m.^7;
                    m = m./max(m(:));
                    val = double(amplitudetodensity(sqrt(S(:,:,1))/2,60,40))/255;
                    HSV(:,:,1) = psi;
                    HSV(:,:,2) = m;
                    HSV(:,:,3) = val;
                    out = hsv2rgb(HSV);
                end
            case 3 % RGB image; nothing to do
            case 4 % Quad-pol
                HH_index=find(strcmpi('H:H',handles.TxRcvPolarizationProc));
                HV_index=find(strcmpi('H:V',handles.TxRcvPolarizationProc));
                VH_index=find(strcmpi('V:H',handles.TxRcvPolarizationProc));
                VV_index=find(strcmpi('V:V',handles.TxRcvPolarizationProc));
                if strcmpi(rep,'AlphaEntropy')
                    tic;
                    [out(:,:,1), out(:,:,2), out(:,:,3)] =ComputeAlphaEntropy(out(:,:,HH_index),...
                        out(:,:,HV_index),out(:,:,VV_index),...
                        out(:,:,VH_index),1);
                    toc;
                    out(:,:,4) = [];
                elseif strcmpi(rep,'Stokes')
                    [S,m,psi] = ComputeStokesQuad(out(:,:,HH_index),out(:,:,HV_index),...
                                out(:,:,VH_index),out(:,:,VV_index),5);
                    psi = psi./max(psi(:));
                    psi = psi.*(175/255);
                    MinM = min(m(:));
                    m = m - MinM;
                    m = m./max(m(:));
                    m = 1 - m;                    
                    val = double(amplitudetodensity(sqrt(S(:,:,1))/2,60,40))/255;
                    HSV(:,:,1) = psi;
                    HSV(:,:,2) = m;
                    HSV(:,:,3) = val;
                    out = hsv2rgb(HSV);
                elseif strcmpi(rep,'Pauli')
                        out=cat(3,abs(out(:,:,HH_index)-out(:,:,VV_index)),...
                            abs(out(:,:,HV_index)+out(:,:,VH_index))/2,...
                            abs(out(:,:,HH_index)+out(:,:,VV_index)));
                end
            
        end
    end
else
    set(handles.image_rep_combo,'Value',find(strcmpi('amplitude',reps)));
end
handles.CurrentImage = out;

if isinteger(out) % Many remap functions won't work on complex ints
    out=single(out);
end

remaps = get(handles.RemapCombo,'String');
remap = remaps{get(handles.RemapCombo,'Value')};

% if nargin>2&&~isempty(remap)
% Apply remap per band
try
    for i=1:size(out,3)
        temp(:,:,i)=feval(remap,out(:,:,i));
    end
catch
    error('ApertureTool:InvalidRemap',['Error applying remap function: ' func2str(remap)]);
end
out=temp; % Use intermediate variable to allow for datatype change
if size(out,3)>1 && isfloat(out)% True color images must be between 0 and 1
    out=out-min(out(isfinite(out)));
    out=out/max(out(:));
end
end

