function [R, RSD, RW] = calculate_linewidth(temp, dims, pps, sdcs, ms)
% Calculates the linewidth of a cell based on the neon and n2 pressures.
%
% pps, sdcs, ms all must be vectors.
% pps -> Partial pressures
% sdcs -> Spin-destruction cross-sections
% ms -> masses
%
% T is in celsius

if(~exist('sdcs', 'var') || sdcs == 0)
   % By default it's neon and nitrogen
   sdcs = [1.9e-23, 1e-22]; % Spin-destruction cross-sections in cm^2
   
   if(length(pps) == 1)
       sdcs = sdcs(1);
   elseif length(pps) ~= 2
       error('Must provide spin-destruction cross sections for more than 2 components.');
   end  
end

if(~exist('ms', 'var') || ms == 0)
    ms = [3.35082*1e-26, 4.65173e-26]; % Mass in kg, Neon, N2
    
    if(length(pps) == 1)
        ms = ms(1);
    elseif(length(pps) ~= 2)
      error('Must provide masses for more than 2 components.');
    end
end

temp = temp + 273.2;

rm = 1.41923e-25;                           % Rubidium mass in kg.
kb = 1.38*1e-23;                            % Boltzmann constant, J/K

pps = pps * 133.3;                          % Pressures should be in pascals, not torr.
n = pps/(kb * 293);                         % Filled at room temperature
n = n * 1e-6;                               % In cm^-3

Mr = (1/rm + 1./ms).^(-1);                  % Reduced masses.
vb = sqrt((8*kb*temp)./(pi*Mr)) * 100;    	% Average atomic velocity in cm/s 

P = 1/2;
q_85 =  (38 + 52 * P^2 + 6 * P^4)/(3 + 10*P^2 + 3 * P^4);
q_87 = (6 + 2*P^2)/(1 + P^2);

RSDB = sum(vb.*sdcs.*n);
RSD85 = (1/q_85)*RSDB;
RSD87 = (1/q_87)*RSDB;

RSD = [RSD85, RSD87];                       % Rate is in s^-1

D = 0;
if(pps(1) > 0)
    D1 = 0.235 * (temp/305).^(3/2) * (760*133.3/pps(1));    % cm^2/s
    D = D1;
end

if(pps(2) > 0)
   D2 = 0.19*(temp/273).^(3/2)*(760*133.3/pps(2));
   if(D ~= 0)
      D = (1./D1 + 1./D2).^(-1); 
   else
       D = D2;
   end
end

RW = D*(sum((pi./dims).^2));

R = RW + RSD;

R = R/(2*pi);




