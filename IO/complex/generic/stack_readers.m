function [ readerobj ] = stack_readers(ro)
%OPEN_STACKED_SET Concatenates file reader objects to output a multiband
% (3-D) array
%
% readerobj = stack_readers(filereaderobjs)
%
% FILEREADEROBJS is a cell array of open_reader file reader objects.  This
% function deals only with reading complex data, so metadata will have to
% be handled separately.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

%% Define object methods
meta=ro{1}.get_meta(); % Assumes all images are the same size
readerobj = chipfun2readerobj(@chipper, [meta.ImageData.NumCols meta.ImageData.NumRows]);
readerobj.close=@close_function;

    function out = chipper(varargin)
        chips=cell(1,length(ro));
        for j=1:length(ro)
            chips{j}=ro{j}.read_chip(varargin{:});
        end
        out = cat(3, chips{:});
    end

    function close_function()
        for j=1:length(ro), ro{j}.close(); end
    end
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////