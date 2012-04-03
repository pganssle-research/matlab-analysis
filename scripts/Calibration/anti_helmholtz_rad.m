function [Br, Bx] = anti_helmholtz_full(x, r, I, n, R)
   % Calculates the radial field as a function of distance from the center  
   % of the coil.
   %
   %  Defaults are I = 1A, n = 4 turns, R = 25.4mm (1")
   %  Units of output are Gauss.
   %
   %  field = anti_helmholtz(x, I, n, R);

   m0 = 4*pi; % Relative permeability in Gauss*mm/Amp;
    
   if(nargin < 2)
       I = 1;      % 1 amp is default
   end
    
   if(nargin < 3)
       n = 4;      % 4 Turns is default
   end
   
   if(nargin < 4)
       R = 25.4;   % In mm - 1" is default
   end
   
   alpha = r/R;
   beta = x/R;
   gamma = x/r;
   Q = ((1+alpha).^2 + beta.^2);
   k = ((4*alpha)./Q).^(1/2);
  
   B0 = I*m0/(2*R);
   
   [K, E] = ellipke(k);
   
   Bx = B0*(1./(pi*Q.^(1/2))).*(E.*((1-alpha.^2-beta.^2)./(Q-4*alpha)) + K);
   Br = B0*(gamma./(pi*Q.^(1/2))).*(E.*((1+alpha.^2+beta.^2)./(Q-4*alpha))-K);
   
end