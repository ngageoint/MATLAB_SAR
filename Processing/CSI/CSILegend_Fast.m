function CSILegend_Fast(filename)
%CSILEGEND_FAST Plot fast time legend
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%get processed bandwidth
reader_obj = open_reader(filename);
if iscell(reader_obj)
    reader_obj = reader_obj{1};
end
meta = reader_obj.get_meta();
reader_obj.close();

if isfield(meta,'CollectionInfo') && isfield(meta.CollectionInfo,'CoreName') 
    IID = meta.CollectionInfo.CoreName(1:min(16,end));
else
    IID = '';
end

if isfield(meta,'Grid') && isfield(meta.Grid,'Row') && ...
        isfield(meta.Grid.Row,'KCtr') && isfield(meta.Grid.Row,'ImpRespBW')
    CenterFreq = meta.Grid.Row.KCtr * SPEED_OF_LIGHT( ) / 2;
    BW = meta.Grid.Row.ImpRespBW * SPEED_OF_LIGHT( ) / 2;
    StartFreq = CenterFreq - BW/2;
    if ~exist('handle','var')
        h = figure;
        set(h,'Position',[100 100 800 200]);
    end
    imagesc(fliplr(repmat(1:256,1000,1)));
    set(gca,'YTick',[]);
    XTick = 0:32:256;
    XTick(1) = 1;
    set(gca,'XTick',XTick);
    Freqs = round(((XTick/256).*BW+StartFreq)./1e6);
    set(gca,'XTickLabel',Freqs);
    xlabel('Frequency (MHz)');
    a = sprintf('Fast-Time CSI Legend for IID: %s',IID);
    title(a);
else
    return;
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////