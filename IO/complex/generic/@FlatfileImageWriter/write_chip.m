function num_bytes_written = write_chip(obj, data, start_indices)
% WRITE_CHIP Writes the given data into the already-opened file.
%
% Note this is a member function of FlatfileImageWriter class and hence
% must be called in the context of an FlatfileImageWriter object.
%
% Inputs:
%       obj:               This is the 'hidden' object-specific handle.  It
%                          is similar to 'this' in Java.
%       data:              The matrix of data to be written
%       start_indices:     1x2 vector ([start_column_index, start_row_index])
%                          giving starting (upper-left-hand) position of
%                          the chip to write in the entire file. (Note:
%                          Definition of "column" and "row" assumes file
%                          is written row-major order.  This is not the
%                          same as MATLAB, which stores arrays in
%                          column-major.)
%
% Outputs:
%       num_bytes_written: The number of bytes written to the file
%                          associated with the object.
%
% Because of the way that MATLAB's fwrite function works, all bytes in the
% data area of the output data file not explictly written to will be filled
% with zeros, even if the file was written to in a random-access fashion.
% There is no need to explicitly clear all values to zero.
%
% Written by: Tom Krauss and Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Error checking
if any(start_indices<1) || any((double(start_indices) + size(data) - 1) > obj.image_size)
    error('FLATFILEIMAGEWRITER_WRITE_CHIP:InvalidSizePosition','Chip cannot be written outside predefined size of image in file.');
end

%% FWRITE doesn't handle complex numbers, so we have to interleave ourselves
if obj.is_complex
    data_interleaved = zeros(size(data).*[2 1],obj.buffer_type);
    data_interleaved(1:2:end,:) = real(data);
    data_interleaved(2:2:end,:) = imag(data);
else
    data_interleaved = data;
end

%% Need to seek to the correct location in the file before writing. The
% conversion of the offsets/sizes to double is to ensure the multiplies
% don't overflow.
element_size = double(obj.data_size * (obj.is_complex + 1)); % Number of bytes for each pixel
row_skip = element_size * (obj.image_size(1) - size(data,1)); % Number of bytes to skip between each row
offset = obj.header_skip + ... % Start of data in file
    (element_size * double(start_indices(2)-1) * double(obj.image_size(1))) + ... % First row
    (element_size * double(start_indices(1)-1)); % First column

%% We might need to write the first line in the file separately.  MATLAB's
% fwrite skips first and THEN writes, so we might not be able to set the
% file pointer to a position where it can skip and then write without
% having a negative position with respect to the beginning of the file.
if (offset - row_skip) < 0
    obj.fseek(obj.FID, offset, 'bof');
    num_bytes_written = obj.data_size * fwrite(obj.FID, data_interleaved(:,1), obj.data_type);
    offset = offset + (element_size * obj.image_size(1)); % Move offset one row
    data_interleaved = data_interleaved(:,2:end); % First row already written
else
    num_bytes_written = 0;
end
obj.fseek(obj.FID, offset - row_skip, 'bof'); % Subtract off one row_skip for skip parameter in fwrite (which skips before writing)

%% Write all data (except maybe first line in file, if written above)
num_bytes_written = num_bytes_written + (obj.data_size * ...
    fwrite(obj.FID, data_interleaved, [num2str(size(data_interleaved,1)) '*' obj.data_type], row_skip));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////