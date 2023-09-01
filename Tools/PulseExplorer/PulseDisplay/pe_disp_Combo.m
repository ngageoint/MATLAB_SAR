function pe_disp_Combo(Parent,PhReader,PulseNum,ChannelNum,CLim,Filename)
%COMBO: Combo RF/Deramp/Deskew or RF/Deskew
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

meta = PhReader.get_meta();
[~, pvp] = PhReader.read_raw(1,1,ChannelNum);
if (pvp.FICRate==0)
    % Original RF signal
    ax(1) = axes('Parent',Parent,'Position',[0.08 0.56 0.85 0.4]);
    pe_disp_RFSignal(ax(1),PhReader,PulseNum,ChannelNum,CLim);
    % Deskewed signal
    ax(2) = axes('Parent',Parent,'Position',[0.08 0.08 0.85 0.4]);
    pe_disp_Deskewed(ax(2),PhReader,PulseNum,ChannelNum,CLim);
    title_string = sprintf('RF/Deskew of Pulse %d: IID: %s',PulseNum,meta.CollectionID.CoreName(1:min(16,end)));
else % Stretch
    % Reramped signal
    ax(1) = axes('Parent',Parent,'Position',[0.08 0.55 0.85 0.43]);
    pe_disp_RFSignal(ax(1),PhReader,PulseNum,ChannelNum,CLim);
    set(ax(1),'XTick',[]);
    xlabel(ax(1),'');
    % Deramped signal
    ax(2) = axes('Parent',Parent,'Position',[0.08 0.34 0.85 0.2]);
    pe_disp_Deramped(ax(2),PhReader,PulseNum,ChannelNum,CLim);
    % Deskewed signal
    ax(3) = axes('Parent',Parent,'Position',[0.08 0.07 0.85 0.2]);
    pe_disp_Deskewed(ax(3),PhReader,PulseNum,ChannelNum,CLim);
    title_string = sprintf('RF/Deramp/Deskew of Pulse %d: IID: %s',PulseNum,meta.CollectionID.CoreName(1:min(16,end)));
end
set(Parent,'Title',title_string);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////