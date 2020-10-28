function phdchip = ApToolFilterPHDChip(phdchip,WindowFilter,FilterRange,FilterAz)   

[ny,nx,nz] = size(phdchip);

for i=1:nz
    switch WindowFilter        
        case 1 %gaussian
            if FilterRange
                win = gausswin(ny);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win,1,nx);
            end
            if FilterAz
                win = gausswin(nx);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win',ny,1);
            end
        case 2 %1/x^4
            if FilterRange
                win = x4win(ny);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win,1,nx);
            end
            if FilterAz
                win = x4win(nx);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win',ny,1);
            end
        case 3 %hamming
            if FilterRange
                win = hamming(ny);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win,1,nx);
            end
            if FilterAz
                win = hamming(nx);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win',ny,1);                
            end
        case 4 %cosine on pedastal
            if FilterRange
                win = cospedwin(ny,0.5);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win,1,nx);
            end
            if FilterAz
                win = cospedwin(nx,0.5);
                phdchip(:,:,i) = phdchip(:,:,i).*repmat(win',ny,1);
            end
    end
end

