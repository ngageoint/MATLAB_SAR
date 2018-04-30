function [HHCof HVCof VHCof VVCof] = ComputePolCoeff(TxAngle,RxAngle,TxEllip,RxEllip)

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

HHCof = complex(cosd(TxAngle)*cosd(TxEllip),-1*sind(TxAngle)*sind(TxEllip))*...
        complex(cosd(RxAngle)*cosd(RxEllip),-1*sind(RxAngle)*sind(RxEllip));
    
HVCof = complex(cosd(TxAngle)*cosd(TxEllip),-1*sind(TxAngle)*sind(TxEllip))*...
        complex(sind(RxAngle)*cosd(RxEllip),cosd(RxAngle)*sind(RxEllip));
    
VHCof = complex(sind(TxAngle)*cosd(TxEllip),cosd(TxAngle)*sind(TxEllip))*...
        complex(cosd(RxAngle)*cosd(RxEllip),-1*sind(RxAngle)*sind(RxEllip));
    
VVCof = complex(sind(TxAngle)*cosd(TxEllip),cosd(TxAngle)*sind(TxEllip))*...
        complex(sind(RxAngle)*cosd(RxEllip),cosd(RxAngle)*sind(RxEllip));   
    
return;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////