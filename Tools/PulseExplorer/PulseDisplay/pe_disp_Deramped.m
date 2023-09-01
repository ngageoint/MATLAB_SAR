function pe_disp_Deramped(Parent,PhReader,PulseNum,ChannelNum,CLim,Filename)
%DERAMPED: Plot of Deramped pulse for PulseExplorer
%
% Author: Tim Cox, NRL, and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

meta = PhReader.get_meta();

[pulse, pvp] = PhReader.read_raw(PulseNum,'All',ChannelNum);
pulse = double(pulse);
pulse(isnan(pulse)) = 0;

sampling_rate =  meta.Channel.Parameters(ChannelNum).Fs;
if (pvp.FICRate==0)
    %need to do fft, then we will create the equivilant of a deramped
    %plot.  This is only to match what stretch data looks like at this
    %stage...we may want to eliminate this option for chirp altogether
    [pulse,delta_t] = reramp(pulse, sampling_rate, - pvp.TxLFM(3)); % Deramp is just "reramp" of the negative chirp rate
    deramped_sampling_rate = 1 / delta_t;
    [ yy_, t_, f_ ] = stft(pulse, deramped_sampling_rate, true);
    toa_extent = numel(pulse)/sampling_rate + pvp.TXmt;
else
    %just stft of pulse
    pulse = vertcat(zeros(1000,1),pulse,zeros(1000,1));
    [ yy_, t_, f_ ] = stft(pulse, sampling_rate, true);        
    toa_extent = abs(sampling_rate/pvp.FICRate);
end

if strcmpi(get(Parent,'type'),'uipanel')
    axes_handle_right = axes('Parent', Parent,...
        'Position',[.08 .07 .85 .9], 'YAxisLocation','right', ...
        'XTick', [], 'Ylim',([-toa_extent toa_extent]/2)*1e3);
    axes_handle = axes('Parent', Parent, ...
        'Position', get(axes_handle_right,'Position'));
    title_string = sprintf('Deramp STFT of Pulse %d: IID: %s',PulseNum,meta.CollectionID.CoreName(1:min(16,end)));
    set(Parent,'Title',title_string,'FontSize',get(axes_handle,'FontSize'));
elseif strcmpi(get(Parent,'type'),'axes')
    axes_handle_right = axes('Parent', get(Parent,'Parent'),...
        'Position',get(Parent,'Position'), 'YAxisLocation','right', ...
        'XTick', [], 'Ylim',([-toa_extent toa_extent]/2)*1e3);
    axes(Parent);
    axes_handle = Parent;
end
imagesc(t_.*1e3, f_./1e6, yy_, 'Parent', axes_handle);
if CLim(2) ~= 0 && CLim(2) ~= 1
    set(axes_handle,'CLim',CLim);
end
ylabel(axes_handle, 'Frequency (MHz)');
ylabel(axes_handle_right,'Time of Arrival (msecs)');
xlabel(axes_handle, 'Sample Time (msec)');
axis(axes_handle,'xy');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////