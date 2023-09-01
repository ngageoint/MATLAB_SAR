function disp_CPHDRFSignal(Parent,PhReader,PulseNum,ChannelNum,CLim,~)
%DESKEWED: Function to display deskewed pulse inside PulseExplorer
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[pulse, pvp] = PhReader.read_cphd(PulseNum,'All',ChannelNum);
pulse = double(pulse);
pulse(isnan(pulse)) = 0;
meta = PhReader.get_meta();

chirp_rate = (1/pvp.aFRR2)*(2/SPEED_OF_LIGHT);
sampling_rate = pvp.SCSS*(numel(pulse)-1);
[yy_, t, f] = stft(deskew_rvp(ifftshift(ifft(fftshift(pulse))), ...
    sampling_rate, chirp_rate, true), sampling_rate, true);
if strcmpi(get(Parent,'type'),'uipanel')
    axes_handle = axes('Parent', Parent, 'Position',[.08 .07 .85 .9]);
    title_string = sprintf('STFT of Pulse %d: IID: %s',PulseNum,meta.CollectionID.CoreName(1:min(16,end)));
    set(Parent,'Title',title_string,'FontSize',get(axes_handle,'FontSize'));
elseif strcmpi(get(Parent,'type'),'axes')
    axes_handle = Parent;
end
imagesc((t-max(t)/2)*1e3,(f+pvp.SC0)/1e9,yy_, 'Parent', axes_handle);
axis xy;
xlabel(axes_handle, 'Sample Time (msec)');
ylabel(axes_handle, 'Frequency (GHz)');
if CLim(2) ~= 0 && CLim(2) ~= 1
    set(axes_handle,'CLim',CLim);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////