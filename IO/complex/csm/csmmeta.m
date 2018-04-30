function csmmeta( filename )
%CSMMETA Saves metadata from a CSM HDF5 file into a temporary text file.
% Utility for writing all metadata in Cosmo SkyMed SCS (single-look complex
% slant) HDF5 file to a text file, which is generally easier to browse
% through than the HDF5, even with the MATLAB HDF5 utilities.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

meta=hdf5info(filename);

fid=fopen('temp.txt','w');
%% Parameters applicable to entire collect
% /
fprintf(fid,'Root\r\n');
for i=1:length(meta.GroupHierarchy.Attributes)
    write_hdf_value(meta.GroupHierarchy.Attributes(i));
end


%% Parameters for each channel
for i=1:length(meta.GroupHierarchy.Groups) % Go through each channel for "ping-pong" mode
    % /S<mm>
    fprintf(fid,'\r\nS%02d\r\n',i);
    for j=1:length(meta.GroupHierarchy.Groups(i).Attributes)
        write_hdf_value(meta.GroupHierarchy.Groups(i).Attributes(j));
    end
    % /S<mm>/SBI
    fprintf(fid,'\r\nS%02d SBI\r\n',i);
    if i==1
        dataset_index=2; % Dataset 2 is complex data (dataset 1 is quicklook)
    else
        dataset_index=1; % Other channels don't have quicklook
    end
    for j=1:length(meta.GroupHierarchy.Groups(i).Datasets(dataset_index).Attributes)
        write_hdf_value(meta.GroupHierarchy.Groups(i).Datasets(dataset_index).Attributes(j));
    end
end

fclose(fid);


function write_hdf_value(attr)
    fprintf(fid,[attr.Shortname ': ']);
    if strcmp(class(attr.Value),'hdf5.h5string')
        fprintf(fid,attr.Value.Data);
    else
        fprintf(fid,num2str(attr.Value(:)'));
    end
    fprintf(fid,'\r\n');
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////