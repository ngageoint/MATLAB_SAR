function [ b ] = polyshift( a, shift )
%POLYSHIFT Compute coefficients of a shifted 1-D polynomial
%
% The input A is a polynomial of length n and represents the polynomial of
% order (n-1):
% a(1) + a(2)*x + a(3)*(x^2) + ... + a(n)*(x^(n-1))
%
% The polynomial shifted by the scalar input SHIFT is:
% a(1) + a(2)*(x+shift) + a(3)*((x+shift)^2) + ... + a(n)*((x+shift)^(n-1))
% =
% b(1) + b(2)*x + b(3)*(x^2) + ... + b(n)*(x^(n-1))
% 
% This function is a direct computation of the coefficients B of the
% shifted polynomial.  "Synthetic division" would be more efficient, but it
% doesn't matter much when only computing a few shifts for small
% polynomials.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

b=zeros(size(a));
for j = 1:numel(a)
    for k=j:numel(a)
        b(j) = b(j) + (a(k)*nchoosek(k-1,j-1)*shift^(k-j));
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////