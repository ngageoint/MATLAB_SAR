function [OverviewName, phdmeta, overviewmeta] = CreateImageOverview(filename,Prefs)
%CREATEIMAGEOVERVIEW Creates an image overview from CPHD data at the
%specified resolution

%generate overview name
[pathstr, fname] = fileparts(filename);
if strcmpi(Prefs.SICDLocation,'Colocated')    
    OverviewName = [pathstr filesep fname '_Overview.nitf'];
else
    OverviewName = [Prefs.SICDLocation filesep fname '_Overview.nitf'];
end

ph_reader = open_ph_reader(filename);
phdmeta = ph_reader.get_meta();
ph_reader.close();

%if specified, use existing Overview file
if Prefs.UseOverview
    if ~isempty(dir(OverviewName))
        reader_obj = open_reader(OverviewName);
        overviewmeta = reader_obj.get_meta();
        reader_obj.close;
        return;
    end
end

%form image at specified resolution
overviewmeta = pfa_file(filename,OverviewName,'Resolution',[Prefs.MaxResolution Prefs.MaxResolution]);

end

