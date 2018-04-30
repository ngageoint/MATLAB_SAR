function [ output_fun ] = sicdweight2fun( grid_rowcol )
%SICDWEIGHT Make a function from a SICD data structure description of a
%complex image weighting
%
%    output_data = sicdweight2fun(grid_rowcol)
%
%       Parameter name    Description
% 
%       grid_rowcol       Either the Grid.Row or Grid.Col SICD field
%                         depending on which direction is being processed.
%                         Should have either the WgtType or WgtFunct
%                         subfields.
%       output_fun        Function handle that generates weighting.  Takes
%                            a single input parameter, which is the number
%                            of elements in the resulting weighting vector.
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

useWgtFunct = false;
% First try to compute function analytically
if isfield(grid_rowcol,'WgtType')
    try % SICD metadata is sometimes misformed
        if ischar(grid_rowcol.WgtType) % True for SICD versions <0.5
            grid_rowcol.WgtType.WindowName = grid_rowcol.WgtType;
        end
        switch upper(grid_rowcol.WgtType.WindowName)
            case 'UNIFORM'
                % We could do this:
                % output_fun = @(x) ones(x,1);
                % Instead we just pass out an empty array as a simple way
                % to let calling functions know that no weighting function
                % was applied.
                output_fun = [];
            case 'HAMMING'
                if ~isfield(grid_rowcol.WgtType,'Parameter') || ...
                        ~isfield(grid_rowcol.WgtType.Parameter,'value')
                    % A Hamming window is defined in many places as a
                    % raised cosine of weight .54, so this is the default.
                    % However some data use a generalized raised cosine and
                    % call it HAMMING, so we allow for both uses.
                    grid_rowcol.WgtType.Parameter.value = num2str(.54);
                end
                output_fun = @(x) raised_cos_fun(x,str2double(grid_rowcol.WgtType.Parameter.value));
            case 'HANNING'
                output_fun = @(x) raised_cos_fun(x,0.5);
            case 'KAISER'
                % We don't use the Mathworks Kaiser function, so we won't be dependent on the Signal Processing Toolbox
                output_fun = @(x) kaiser_nosptb(x,str2double(grid_rowcol.WgtType.Parameter.value));
            case 'TAYLOR'
                params = grid_rowcol.WgtType.Parameter; % To shorten code
                nbar = str2double(params{cellfun(@(x) strcmpi(x.name,'NBAR'),params)}.value);
                sll = str2double(params{cellfun(@(x) strcmpi(x.name,'SLL'),params)}.value);
                sll = -abs(sll); % A couple conventions for how this may be populated, but only one sign makes sense for MATLAB function
                output_fun = @(x) taylorwin(x,nbar,sll)/max(taylorwin(x,nbar,sll));
            otherwise
                useWgtFunct = true;
        end
        if ~isempty(output_fun)
            % Run once just to make sure the function we created doesn't throw error
            feval(output_fun,2);
        end
    catch
        useWgtFunct = true;
    end
else
    useWgtFunct = true;
end
% If analytic approach didn't work, use sampled data
if useWgtFunct
    if ~isfield(grid_rowcol,'WgtFunct')
        error('SICDWEIGHT2FUN:INSUFFICIENT_METADATA','Insufficient metadata to determine weighting function.');
        % Alternative for calling functions, if they catch this error, is
        % to determine weighting from the complex data itself.
    else
        output_fun = @(x) interpft(grid_rowcol.WgtFunct,x);
    end
end

end

% Reproduces hamming functionality from the MATLAB Signal Processing
% Toolbox, but allows for arbitrary coefficients of raised cosine.
function w = raised_cos_fun(n, coef)
    w = coef - (1-coef)*cos(2*pi*(0:ceil(n/2)-1)'/(n-1));
    if ~rem(n,2)
        w = [w; w(end:-1:1)];
    else
        w = [w; w(end-1:-1:1)];
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////