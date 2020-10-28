function ApToolupdatePHD_axesLabels(handles)
%update PHD axes labels

xlim = get(handles.phd,'xlim');
ylim = get(handles.phd,'ylim');

%setup axes labels
if isfield(handles,'phdXaxis')
    dxlim = diff(xlim);
    xticks = round(xlim(1))+[0 ceil(dxlim/4-1) ceil(dxlim/2-1) ...
        ceil(dxlim/4*3-1) dxlim];
    xticks(xticks < xlim(1)) = ceil(xlim(1));
    xticks(xticks > xlim(2)) = floor(xlim(2));
    xticklabels = cellfun(@(x) sprintf('%0.2g',x),num2cell(handles.phdXaxis(xticks)),'uniformoutput',false);
    xticklabels{3} = '0'; % Should be very close to zero, so this makes display cleaner.
    set(handles.phd,'xtick',xticks,'xticklabel',xticklabels);
    xlabel(handles.phd,'Polar Angle (degrees)','fontweight','bold')
else
    set(handles.phd, 'xtick', []);
end
if isfield(handles,'phdYaxis')
    dylim = diff(ylim);
    yticks = round(ylim(1))+[0 ceil(dylim/4-1) ceil(dylim/2-1) ...
        ceil(dylim/4*3-1) dylim];
    yticks(yticks < ylim(1)) = ceil(ylim(1));
    yticks(yticks > ylim(2)) = floor(ylim(2));
    set(handles.phd,'ytick',yticks,'yticklabel',...
        cellfun(@(x) sprintf('%0.4g',x),num2cell(handles.phdYaxis(yticks)),'uniformoutput',false));
    ylabel(handles.phd,sprintf('Frequency (%s)',handles.freqUnits),'fontweight','bold')
else
    set(handles.phd, 'ytick', []);
end

