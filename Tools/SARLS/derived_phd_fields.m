function meta = derived_phd_fields(meta,nbdata)
%DERIVED_PHD_FIELDS - Function fills in SICD meta fields not populated for
%phd data. These fields aren't explicitly known until image formation, but
%with some assumptions (i.e. all pulses processed, full res), we can
%populate these fields for display in tools.  It's often nice to know what
%these fields will generally be if an image is formed.

%set up fields to call derived_sicd_fields to populate SCPCOA
if isfield(nbdata,'SRPPos')
    meta.GeoData.SCP.ECF.X = nbdata.SRPPos(2,1);
    meta.GeoData.SCP.ECF.Y = nbdata.SRPPos(2,2);
    meta.GeoData.SCP.ECF.Z = nbdata.SRPPos(2,3);
elseif isfield(meta,'ReferenceGeometry') && isfield(meta.ReferenceGeometry,'CRP')
    meta.GeoData.SCP = meta.ReferenceGeometry.CRP;
end

meta.SCPCOA.ARPPos.X = nbdata.RcvPos(2,1);
meta.SCPCOA.ARPPos.Y = nbdata.RcvPos(2,2);
meta.SCPCOA.ARPPos.Z = nbdata.RcvPos(2,3);

meta.SCPCOA.ARPVel.X = (nbdata.RcvPos(3,1)-nbdata.RcvPos(2,1))/(nbdata.RcvTime(3)-nbdata.RcvTime(2));
meta.SCPCOA.ARPVel.Y = (nbdata.RcvPos(3,2)-nbdata.RcvPos(2,2))/(nbdata.RcvTime(3)-nbdata.RcvTime(2));
meta.SCPCOA.ARPVel.Z = (nbdata.RcvPos(3,3)-nbdata.RcvPos(2,3))/(nbdata.RcvTime(3)-nbdata.RcvTime(2));

meta.Timeline.CollectDuration = nbdata.TxTime(end) - nbdata.TxTime(1);

%fill in Timeline IPP for TStart, End and Pulses.  This will allow for PRF
%computation
meta.Timeline.IPP.Set.TStart = nbdata.TxTime(1);
meta.Timeline.IPP.Set.TEnd = nbdata.TxTime(end);
meta.Timeline.IPP.Set.IPPStart = 1;
meta.Timeline.IPP.Set.IPPEnd = meta.Data.Channel.NumVectors;

meta = derived_sicd_fields(meta);

end

