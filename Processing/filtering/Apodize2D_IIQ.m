%--------------------------------------------------
% S = Apodize2D_IIQ(C, k_h, k_v)
%--------------------------------------------------
%
% Returns a complex image which is created by applying Spatially 
% Variant Apodization (SVA) to the supplied input complex image.
%
% 2D SVA 
%   - Two dimensions simultaneously
%   - I,Q separately
%   - Two dimensions uncoupled
%   
% see equation (41) (plus conditions)
%
% Inputs:
%       C:     The complex image to which SVA is to be applied
%       k_h:   The image sample rate in the horizontal direction
%       k_v:   The image sample rate in the vertical direction
%
% Note: This is in no way optimized (for speed or memory usage).  Rather it
% is meant to make the implementation as obvious as possible.  As such, it
% attempts to map the algorithms, equations, and variables as used in the
% paper "Nonlinear Apodization for Sidelobe Control in SAR Imagery",
% Stankwitz et al.

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

function im_s = Apodize2D_IIQ(im, k_h, k_v)
    yi = DoRealApodize(real(im),k_h,k_v);
    yq = DoRealApodize(imag(im),k_h,k_v);

    im_s = yi + j*yq;
end

function g_prime = DoRealApodize(g, k_h, k_v)
   [rows,cols] = size(g);
   g_prime = zeros(rows,cols);
   for m=1+k_h:rows-k_h
      for n=1+k_h:cols-k_v
          Q_m = g(m-k_h,n) + g(m+k_h,n);
          Q_n = g(m,n-k_v) + g(m,n+k_v);
          P = g(m-k_h,n-k_v) + ...
              g(m+k_h,n+k_v) + ...
              g(m-k_h,n+k_v) + ...
              g(m+k_h,n-k_v);

             w_m = 0.0;
             w_n = 0.5;
             g_p1 = g(m,n) + w_m*w_n*P + w_m*Q_m + w_n*Q_n;

             w_m = 0.5;
             w_n = 0.0;
             g_p2 = g(m,n) + w_m*w_n*P + w_m*Q_m + w_n*Q_n;

             w_m = 0.5;
             w_n = 0.5;
             g_p3 = g(m,n) + w_m*w_n*P + w_m*Q_m + w_n*Q_n;

%          if (sign(g_p1) ~= sign(g(m,n)) | ...
%              sign(g_p2) ~= sign(g(m,n)) | ...
%              sign(g_p3) ~= sign(g(m,n)))
          if (g_p1*g(m,n)<0 || ...
              g_p2*g(m,n)<0 || ...
              g_p3*g(m,n)<0)
             g_prime(m,n) = 0.0;
          else
             data = [g(m,n), g_p1, g_p2, g_p3];
             [d_v,d_i] = min(abs(data));
             %g_prime(m,n) = min( [g(m,n), g_p1, g_p2, g_p3] );
             g_prime(m,n) = data(d_i);
          end
      end
   end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////