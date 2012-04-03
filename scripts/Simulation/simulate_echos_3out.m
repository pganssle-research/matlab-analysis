function out = simulate_echos_3out(time, x_offset, y_offset, z_calib)
    % Simulates  the effect of a quadrature pulse with z readout with imperfect conditions.
    % Zero relaxation, using hydrogens, readout along z.
    % 
    % Time = time between the echo pulses, in ms. Default is 200
    %
    % x, y and z offsets are amount that the pulses are wrong in radians.
    %
    % Default values are: 0, 0, 0.53G
    %
    % out = 2D vector reads out x, y and z components at each step.
    % Step 1: After the first x pi/2.
    % Step 2: After the first evolution period under z.
    % Step 3: After the x pi
    % Step 4: After the y pi
    % Step 5: After the x pi
    % Step 6: After the final evolution under z
    % Step 7: After the final x pi/2.
    %
    % Usage: out = simulate_echos_3out(time, x_offset, y_offset, z_calib)
   

    if nargin == 0
       time = 200;
    end

    if nargin < 2
        x_offset = 0;
    end
    
    if nargin < 3
        y_offset = 0;
    end
      
    if nargin < 4
        z_calib = 0.53;
    end
    
    ang = 26.75 * z_calib *time; % Angle swept out in "time" 26.75 kHz/Gauss * Gauss * time in ms
    
    % Set up the x y and z coordinates.
    x = [1; 0; 0];
    y = [0; 1; 0];
    z = [0; 0; 1];
    
    pi_x = pi+x_offset;
    pi2_x = pi/2 + x_offset/2;
    
    pi_y = pi+y_offset;
    pi2_y = pi/2 + y_offset/2;
    
    % Set up the vector:
    v = z;
    
    % First pulse, pi/2 x to put it in the xy plane.
    v = RotVecArAxe(v, x, pi2_x);
    out(1, :) = v;
    
    % Then let it evolve for some time, record.
    v = RotVecArAxe(v, z, ang);
    out(2, :) = v; 
   
    % Immediate pi pulse along y, record.
    
    v = RotVecArAxe(v, x, pi2_x);  
    v = RotVecArAxe(v, y, pi2_y);
    v = RotVecArAxe(v, x, pi2_x);
    v = RotVecArAxe(v, y, pi2_y);
    v = RotVecArAxe(v, x, pi2_x);  
    
       
    % Evolve again for a bit, record.
    v = RotVecArAxe(v, z, ang);
    out(4, :) = v;
    
    % Finally, the last pi/2 pulse.
    v = RotVecArAxe(v, x, pi2_x);
    out(5, :) = v;
     
    
    
function B = RotVecArAxe(A,L,Phi)
    % Function B = RotVecArAxe(L,Phi,A) rotate 
    % the vector A(3,1) around the axe L(3,1) 
    % into angle Phi radian counterclockwise.
    % Author: Sergiy Iglin
    % e-mail: iglin@kpi.kharkov.ua
    % or: siglin@yandex.ru
    % personal page: http://iglin.exponenta.ru

    %========== Checking of datas ==================
    strerr='Is needed 3 input parameters: vector A(3,1), axe L(3,1) and angle Phi!';
    if nargin<3,
      error(strerr);
    end
    if ~isnumeric(A)
      error(strerr);
    else
      A=A(:);
      if length(A)<3,
        error(strerr);
      end
      A=A(1:3);
    end
    if ~isnumeric(L)
      error(strerr);
    else
      L=L(:);
      if length(L)<3,
        error(strerr);
      end
      L=L(1:3);
      L0=L/norm(L); % orth
    end
    if ~isnumeric(Phi)
      error(strerr);
    else
      Phi=Phi(1);
    end

    %============= Solution ================
    cphi=cos(Phi);
    B=A*cphi+(A'*L0)*(1-cphi)*L0+cross(L0,A)*sin(Phi);
    return
