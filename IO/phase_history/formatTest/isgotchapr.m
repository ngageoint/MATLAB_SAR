function boolout = isgotchapr( filename )
%ISGOTCHAPR Checks to see whether FILENAME is a file in a directory that
%matches the format of the AFRL 2D/3D SAR Volumetric public release dataset
%
% The GOTCHA public release volumetric data is available by request at:
% https://www.sdms.afrl.af.mil/index.php?collection=gotcha
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

[path, name, ext] = fileparts(filename);
boolout = regexp([name ext],'data_3dsar_pass[1-9]_az[0-9]{3}_[HV]{2}.mat','once');
% This only works for MATLAB 2011b
% if boolout
%     try % Check that file is a valid MAT file with required data
%         info = whos(matfile(filename));
%         boolout = (length(info)==1)&&strcmp(info.name,'data')&&strcmp(info.class,'struct');
%     catch
%         boolout = false;
%     end
% end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////