function im = UpsampleImage(img, line_factor, elem_factor)
% Returns an FFT upsampled image.
%
% Returns an image which is upsampled by the inverse FFTing, zero padding
% by the supplied amount (line_factor and elem_factor) and FFTing back.
%
% Usage:
%     im = UpsampleImage(img, line_factor, elem_factor)
%
% Example:
%     im = UpsampleImage(some_image data, 1.5, 1.5);
%
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

    [rows,cols] = size(img);

    % Inverse FFT
    sp = fftshift(ifft2(fftshift(img))) / (rows*cols);
    %%%sp = ifft2(img) / (rows*cols);

    padded_data = zeros( floor(rows*line_factor), floor(cols*elem_factor) );
    if (line_factor < 1  &&  elem_factor < 1)
        row_start  = floor((rows-(rows*line_factor))/2);
        row_length = size(padded_data,1);
        col_start  = floor((cols-(cols*elem_factor))/2);
        col_length = size(padded_data,2);
        padded_data = sp(row_start:(row_start+row_length), ...
                         col_start:(col_start+col_length));
    else
        % Zero pad around the edges
        off_line = floor((rows*line_factor - rows)/2.0);
        off_elem = floor((cols*elem_factor - cols)/2.0) ;
        padded_data(off_line:off_line+rows-1, off_elem:off_elem+cols-1) = sp;
    end

    % Then just FFT back to get the image.
    im = fftshift(fft2(fftshift(padded_data))) * (rows*cols);
    %%%im = fft2(padded_data) / (rows*cols);
end
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
