function [wbvectors, nbdata, meta] = read_cphd_data(filename, varargin)
% READ_CPHD_DATA Generic function to read compensated phase history data
% from ANY file format of phase history data (not just the file format
% known as CPHD)
%
% See documentation in OPEN_PH_READER for syntax on how to call this
% function with input and output arguments.  Syntax for READ_CPHD_DATA is
% the same as the read_cphd method produced by OPEN_PH_READER, only with a
% filename as the first argument.
%
% This is really just a wrapper to the phase history reader objects created
% by OPEN_PH_READER.  This allows users to avoid the objected-oriented
% syntax of those readers, but still provides sensor- and format-
% independent reading of phase history data.  However, this function
% requires that for every call to it, a file be opened, its metadata
% completely parsed, and the file be closed again for each call.
% Therefore, if one will be repeatedly reading pulse data from the same
% file, it would be much faster to use the OPEN_PH_READER function.  For
% single reads (or only occasional reads) though, READ_CPHD_DATA should
% work just fine.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

file_reader_object=open_ph_reader(filename);
try % File must still be closed, even if read fails
    [wbvectors, nbdata]=file_reader_object.read_cphd(varargin{:});
    meta=file_reader_object.get_meta();
catch
    file_reader_object.close();
    rethrow(lasterror);
end
file_reader_object.close();

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////