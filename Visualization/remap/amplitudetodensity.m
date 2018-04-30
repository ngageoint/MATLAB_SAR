function D = amplitudetodensity(A, Dmin, Mmult, data_mean)
% AMPLITUDETODENSITY Convert to density data for remap
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

  eps = 1E-5;
  
  if nargin<2
      Dmin = 30;
  end
  if nargin<3
      Mmult = 40;
  end
  
  A=abs(single(A)); % Most of this basic math won't work in MATLAB on complex integers
  if nargin<4
      data_mean=mean(A(isfinite(A(:))));
  end
  Cl = 0.8*data_mean;
  Ch = Mmult*Cl;
  m = (255-Dmin)/log10(Ch/Cl);
  b = Dmin - m*log10(Cl);
  
  D = m*log10(max(A,eps)) + b;
end
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

