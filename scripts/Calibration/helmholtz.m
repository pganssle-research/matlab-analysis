function field = helmholtz(x, I, n, R)
   % Calculates the field as a function of distance from the origin in an
   % helmholtz coil.
   %
   %  Defaults are I = 1A, n = 8 turns, R = 36.5125mm (1.4375")
   %  Units of output are Gauss.
   %
   %  field = helmholtz(x, I, n, R);

   m0 = 4*pi; % Relative permeability in Gauss*mm/Amp;
    
   if(nargin < 2)
       I = 1;      % 1 amp is default
   end
    
   if(nargin < 3)
       n = 8;      % 4 Turns is default
   end
    
   if(nargin < 4)
       R = 36.5125;   % In mm - 1" is default
   end
    
   field = (n*m0*I*(R^2)/2)*((R^2 + ((R/2) + x).^2).^(-3/2) + (R^2 + ((R/2) - x).^2).^(-3/2));
    
end