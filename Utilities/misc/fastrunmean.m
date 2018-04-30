function [outmat] = fastrunmean(inmat,win,ptype);
%--------------------------------------------------------------------------
% [outmat] = fastrunmean(inmat,win,ptype);
%
% computes the running mean of 'inmat' over a window of size 
% 'win'.  The ptype is the padding type, which controls how the edges are 
% padded.  An efficient recursive algorithm is used that is independent of the
% window size so that the cost is only proportional to the size of the
% input matrix.
%
% Input
% inmat   - a 2D or 3D matrix 
% win     - a matix of window lengths in each dimension.  Windows should be
%           odd in size!
% ptype   - a string decribing how to treat edge padding
%           'zeros' - pad with zeros
%           'mean'  - pad with the mean of inmat
%           note: plan to add more options when I have time
%
% Output
% outmat  - boxcar filtered inmat.
%
% Example and speed test
%
% m2D = randn([100 100])
% m3D = randn([100 100 100])
%
% 2D
% d = fastrunmean(m2D,[51 51],'zeros');
% 3D
% d = fastrunmean(m3D,[51 21 73],'zeros');
%
% for n=1:100
% tic;d=fastrunmean(m2D,[51 51],'zeros');t(n)=toc;
% end
% fprinft('Average time over 100 calls = %d',mean(t))
%
% for n=1:100
% tic;d=fastrunmean(m3D,[51 51 51],'zeros');t(n)=toc;
% end
% fprinft('Average time over 100 calls = %d',mean(t))
%
% You can alter the size of the windows to see that t is independent of the 
% window size.  The time is averaged to allow for flucuations in computer
% performace.  The average time will give you a feel for the speed.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Author  : Neil Hodgson
% Date    : 10 April 2008 15:00
% Contact : hodgson.neil@yahoo.co.uk
% -------------------------------------------------------------------------

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% -------------------------------------------------------------------------
% Error checking
% -------------------------------------------------------------------------

if ndims(inmat) ~= length(win)
    error('The number of entries in win does match ndims(inmat)')
end

evens = find(rem(win/2,2)==0);
if isempty(evens)==0
    error('One or more the window sizes are even.') 
end
        
% -------------------------------------------------------------------------

if length(win)==2
    
    [nx ny] = size(inmat);
    
    nxy     = nx*ny;
    winx    = win(1);
    winy    = win(2);
    
    hwinx = (winx - 1)/2;
    hwiny = (winy - 1)/2;
    
    % pad with zeros
    if strcmp(ptype,'zeros')
        
        inmatpad = zeros([nx+(2*hwinx) ny+(2*hwiny)]);
        
    % pad with the mean   
    elseif strcmp(ptype,'mean')
    
        inmatpad = mean(mean(inmat)).*...
            ones([ nx+(2*hwinx) ny+(2*hwiny)]);
    end
    
    a = (hwinx+1):(nx+hwinx);
    b = (hwiny+1):(ny+hwiny);
    
    inmatpad(a,b) = inmat;
    
    cx = zeros([nx ny+winy-1]);
    cy = zeros([nx ny]);
    
    cx(1,:) = sum(inmatpad(1:winx,:),1);
    
    for n=2:nx
        cx(n,:) = cx(n-1,:) - inmatpad(n-1,:) + inmatpad(n+winx-1,:);
    end

    cy(:,1) = sum(cx(:,1:winy),2);
    
    for n=2:ny
        cy(:,n) = cy(:,n-1) - cx(:,n-1) + cx(:,n+winy-1);
    end
 
    outmat = cy./(winx*winy);

elseif length(win)==3
    
    [nx ny nz] = size(inmat);
    
    nxyz    = nx*ny*nz;
    winx    = win(1);
    winy    = win(2);
    winz    = win(3);
    
    hwinx = (winx - 1)/2;
    hwiny = (winy - 1)/2;
    hwinz = (winz - 1)/2;
    
    % pad with zeros first test
    if strcmp(ptype,'zeros')
        
        inmatpad = zeros([nx+(2*hwinx) ny+(2*hwiny) nz+(2*hwinz)]);
        
    elseif strcmp(ptype,'mean')
    
        inmatpad = mean(mean(mean(inmat))).*...
            ones([ nx+(2*hwinx) ny+(2*hwiny) nz+(2*hwinz)]);
    end
        
    a = (hwinx+1):(nx+hwinx);
    b = (hwiny+1):(ny+hwiny);
    c = (hwinz+1):(nz+hwinz);
        
    inmatpad(a,b,c) = inmat;
    
    cx = zeros([nx ny+winy-1 nz+winz-1]);
    cy = zeros([nx ny nz+winz-1]);
    cz = zeros([nx ny nz]);
    
    cx(1,:,:) = sum(inmatpad(1:winx,:,:),1);
    
    for n=2:nx
        cx(n,:,:) = cx(n-1,:,:) - inmatpad(n-1,:,:) + inmatpad(n+winx-1,:,:);
    end
    
    cy(:,1,:) = sum(cx(:,1:winy,:),2);
    
    for n=2:ny
        cy(:,n,:) = cy(:,n-1,:) - cx(:,n-1,:) + cx(:,n+winy-1,:);
    end

    cz(:,:,1) = sum(cy(:,:,1:winz),3);
    
    for n=2:nz
        cz(:,:,n) = cz(:,:,n-1) - cy(:,:,n-1) + cy(:,:,n+winz-1);
    end
    
    outmat = cz./(winx*winy*winx);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////