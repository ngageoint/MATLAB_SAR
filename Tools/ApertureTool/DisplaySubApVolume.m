function DisplaySubApVolume(phd,ZpLimsAz,ZpLimsRn,NumFrames,ApFrac,SlowFlag,WindowFilter,Bounds,settings)
%DISPLAYSUBAPVOLUME Computes and displays a 3D volume of the subap data

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////


[ny,nx] = size(phd);

ZpWidthAz = diff(ZpLimsAz);
ZpWidthRn = diff(ZpLimsRn);
ApWidth = round(ZpWidthAz*ApFrac);
ApHeight = round(ZpWidthRn*ApFrac);

DataVolume = zeros(NumFrames,nx,ny);

for CurFrame=1:NumFrames   
    if SlowFlag    
        xmin = max(1,ZpLimsAz(1) + ...
            round((CurFrame-1)*(ZpWidthAz-ApWidth)/(NumFrames-1)));
        xmax = xmin+ApWidth;
        Ap = zeros(ny,nx);
        phdchip = phd(:,xmin:xmax);
        [ny2,nx2] = size(phdchip);
        switch WindowFilter
            case 1 %gaussian
                win = gausswin(nx2);
                phdchip = phdchip.*repmat(win',ny,1);
            case 2 %1/x^4
                win = x4win(nx2);
                phdchip = phdchip.*repmat(win',ny,1);
            case 3 %hamming
                win = hamming(nx2);
                phdchip = phdchip.*repmat(win',ny,1);   
            case 4 %cosine on pedastal
                win = cospedwin(nx2,0.5);
                phdchip = phdchip.*repmat(win',ny,1);
        end
        Ap(:,xmin:xmax) = phdchip;
        DataVolume(CurFrame,:,:) = flipud(abs(fft2(Ap')));
    else
        ymin = max(1,ZpLimsRn(1) + ...
            round((CurFrame-1)*(ZpWidthRn-ApHeight)/(NumFrames-1)));
        ymax = ymin+ApHeight;
        Ap = zeros(ny,nx);
        phdchip = phd(ymin:ymax,:);
        [ny2,nx2] = size(phdchip);
        switch WindowFilter
            case 1 %gaussian
                win = gausswin(ny2)';
                phdchip = phdchip.*repmat(win',1,nx);
            case 2 %1/x^4
                win = x4win(ny2)';
                phdchip = phdchip.*repmat(win',1,nx);
            case 3 %hamming
                win = hamming(ny2)';
                phdchip = phdchip.*repmat(win',1,nx);   
            case 4 %cosine on pedastal
                win = cospedwin(ny2,0.5)';
                phdchip = phdchip.*repmat(win',1,nx);
        end
        Ap(ymin:ymax,:) = phdchip;
        DataVolume(NumFrames-CurFrame+1,:,:) = flipud(abs(fft2(Ap')));
    end

end

voldata = DataVolume;
cutoff = prctile(voldata(:),settings.Percentile);
voldata(voldata<cutoff) = 0;
voldata = voldata - cutoff;

voldata = double(densityremap(voldata));

figure;
vol3d('CData',voldata,'texture','3D','YData',Bounds);
alphamap('rampup');
alphamap(settings.AlphaFactor .* alphamap);

set(gca,'View',[30 30]);
xlabel('Cross Range Position');
if SlowFlag
    ylim(Bounds);
    ylabel('Time(sec)');
else
    ylabel('Frequency(MHz)');
end

zlabel('Range Position');
xlim([0 nx]);
zlim([0 ny]);

settings.CMapCutoff = round(settings.CMapCutoff);
cmap = jet(256);
cmap(1:settings.CMapCutoff,:) = repmat([1 1 1],settings.CMapCutoff,1);
cmap(settings.CMapCutoff:end,:) = jet(256-settings.CMapCutoff+1);
colormap(cmap);

