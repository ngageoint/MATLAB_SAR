function pe_disp_Deskewed(Parent,PhReader,PulseNum,ChannelNum,CLim,Filename)
%DESKEWED: Function to display deskewed pulse inside PulseExplorer
%
% Author: Tim Cox, NRL, and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[pulse, pvp] = PhReader.read_raw(PulseNum,'All',ChannelNum);
pulse = double(pulse);
pulse(isnan(pulse)) = 0;
meta = PhReader.get_meta();

sampling_rate =  meta.Channel.Parameters(ChannelNum).Fs;
if (pvp.FICRate==0)
    % This code assumes that the chirp that passes through DC at the
    % center of the receive window is the reference time-of-arrival. In
    % processing, this reference would be the return corresponding to
    % the motion comp point, but just for the purpose of visualization,
    % this works fine.
    pulse_time_ref = ifftshift(pulse); % Put center of window at time zero
    pulse_fft = fft(pulse_time_ref)/sqrt(numel(pulse_time_ref)); % Normalize so energy is the same
    deskewed_pulse = reramp(fftshift(pulse_fft), ...
        size(pulse,1)/sampling_rate, 1/pvp.TxLFM(3));
    deskewed_sampling_rate = sampling_rate*size(deskewed_pulse,1)/size(pulse,1);
    yy_ = stft(deskewed_pulse, deskewed_sampling_rate, true);
    fc = meta.Channel.Parameters(ChannelNum).F0Ref + pvp.DFIC0;
    freq_axis = linspace(fc-sampling_rate/2,fc+sampling_rate/2,size(yy_,2));
    % Note that the time-of-arrival extent computed here is not the
    % range of "useful" time-of-arrivals, that is, full bandwidth pulse
    % returns (This would be wfp.RcvWindowLength - wfp.TxPulseLength).
    % The toa_extent computed here also includes the partial bandwidth
    % returns out of the deskew, which is useful for display, but not
    % necessarily for processing.
    toa_extent = numel(pulse)/sampling_rate + pvp.TXmt;
else
    [deskewed_pulse, t_start] = deskew_rvp(pulse, sampling_rate, pvp.FICRate, true);   
    yy_ = stft(deskewed_pulse, sampling_rate, true);   
    freq_axis = pvp.RefFreq + pvp.DFIC0 + ...
        ((0:length(deskewed_pulse))' * pvp.FICRate / sampling_rate) + ...
        (t_start * pvp.FICRate);
    toa_extent = sampling_rate/pvp.FICRate;
end
toa_axis = -((1:size(yy_,1))-ceil(size(yy_,1)/2))*toa_extent/size(yy_,1);

if strcmpi(get(Parent,'type'),'uipanel')
    axes_handle = axes('Parent', Parent, 'Position',[.08 .07 .85 .9]);
    title_string = sprintf('Deskew STFT of Pulse %d: IID: %s',PulseNum,meta.CollectionID.CoreName(1:min(16,end)));
    set(Parent,'Title',title_string,'FontSize',get(axes_handle,'FontSize'));
elseif strcmpi(get(Parent,'type'),'axes')
    axes_handle = Parent;
end
imagesc(freq_axis./1e9, toa_axis.*1e3, yy_, 'Parent', axes_handle);
ylabel(axes_handle, 'Time of Arrival (msec)');
xlabel(axes_handle, 'Frequency (GHz)');
if CLim(2) ~= 0 && CLim(2) ~= 1
    set(axes_handle,'CLim',CLim);
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////