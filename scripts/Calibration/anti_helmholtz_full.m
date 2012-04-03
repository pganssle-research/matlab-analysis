function [Br, Bx] = anti_helmholtz_full(x, r, I, n, R)
   % Calculates the radial field as a function of distance from the center  
   % of the coil.
   %
   %  Defaults are I = 1A, n = 4 turns, R = 25.4mm (1")
   %  Units of output are Gauss.
   %
   %  field = anti_helmholtz(x, I, n, R);

   m0 = 4*pi; % Relative permeability in Gauss*mm/Amp;
    
   if(nargin < 3)
       I = 1;      % 1 amp is default
   end
    
   if(nargin < 4)
       n = 8;      % 4 Turns is default
   end
   
   if(nargin < 5)
       R = 25.4;   % In mm - 1" is default
   end
   
   % Ensure that these are row vectors
   sx = size(x);
   sr = size(r);
   
   if(sx(1) > 1)
       x = x';
   end
    
   if(sr(1) > 1)
       r = r';
   end  
   
   % Turn these into matrices of the same size
   x = x'*ones(sr);
   r = ones(sx)'*r;
   
   % Now we do the actual calculations
   alpha = r/R;
   beta = x/R;
   gamma = x./r;
   Q = ((1+alpha).^2 + beta.^2);
   k = ((4*alpha)./Q).^(1/2);
  
   B0 = I*m0/(2*R);
   
   [K, E] = ellipke(k);
   
   Bx = B0*(1./(pi*Q.^(1/2))).*(E.*((1-alpha.^2-beta.^2)./(Q-4*alpha)) + K);
   Br = B0*(gamma./(pi*Q.^(1/2))).*(E.*((1+alpha.^2+beta.^2)./(Q-4*alpha))-K);
   
end