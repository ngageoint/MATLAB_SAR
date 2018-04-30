function structout = setstructfields(structin, newfields)
%SETSTRUCTFIELDS Set fields of a structure using another structure
%   SETSTRUCTFIELDS(STRUCTIN, NEWFIELDS) Set fields of STRUCTIN using
%   another structure NEWFIELDS fields.  If fields exist in STRUCTIN
%   but not in NEWFIELDS, they will not be changed.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.3.4.3 $  $Date: 2007/12/14 15:16:12 $
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%narginchk(2,2);

% Set the output to the input
structout = structin;

% If newfields is not a structure return
if ~isstruct(newfields), return; end

% Loop over all the fields of newfields and assign them to the output
fields = fieldnames(newfields);
for i = 1:length(fields),
    value = newfields.(fields{i});
    if isstruct(value) && isfield(structout, fields{i}) && isstruct(structout.(fields{i}))
        structout.(fields{i}) = setstructfields(structout.(fields{i}), value);
    else
        structout.(fields{i}) = newfields.(fields{i});
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////