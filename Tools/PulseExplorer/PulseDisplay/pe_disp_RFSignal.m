function pe_disp_RFSignal(Parent,PhReader,PulseNum,ChannelNum,CLim,Filename)
%RFSIGNAL: Function to display RF signal (reramped pulse for strech mode,
%raw for chirp) inside PulseExplorer
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
    [ yy_, t_, f_ ] = stft(pulse, sampling_rate, true);
    f_ = f_ + meta.Channel.Parameters(ChannelNum).F0Ref + pvp.DFIC0 - (sampling_rate/2); % Make real frequency
else       
    [cc, dt] = reramp(pulse, sampling_rate, pvp.FICRate);
    reramped_sample_rate = 1/dt;        
    [ yy_, t_, f_ ] = stft(cc, reramped_sample_rate, true);
    f_ = f_ + ... % Add offset to make real receive frequency
        meta.Channel.Parameters(ChannelNum).F0Ref + pvp.DFIC0 + ...
        (pvp.FICRate*size(pulse,1)/(sampling_rate*2))-... % Deramp frequency at center of receive window (which reramp.m places at DC)
        (reramped_sample_rate/2); % Center frequency in reramped bandwidth (as it comes out of stft.m)
end

%resize image so that it is square
if strcmpi(get(Parent,'type'),'uipanel')
    axes_handle = axes('Parent',Parent,'Position',[.08 .07 .85 .9]); 
    title_string = sprintf('RF STFT of Pulse %d: IID: %s',PulseNum,meta.CollectionID.CoreName(1:min(16,end)));
    set(Parent,'Title',title_string,'FontSize',get(axes_handle,'FontSize'));
elseif strcmpi(get(Parent,'type'),'axes')
    axes_handle = Parent;
end
pos = getpixelposition(axes_handle);

if exist('imresize') == 2 %#ok<*EXIST>
    yy_ = imresize(yy_, pos([4 3]));
end

t_ = linspace(t_(1),t_(end),size(yy_,2));

imagesc(t_.*1e3, f_./1e9, yy_, 'Parent', axes_handle);
if CLim(2) ~= 0 && CLim(2) ~= 1
    set(axes_handle,'CLim',CLim);
end
ylabel(axes_handle, 'Frequency (GHz)');
xlabel(axes_handle, 'Sample Time (msec)');
axis(axes_handle,'xy');

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////