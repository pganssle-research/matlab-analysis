function [out, t] = simulate_pulses_3out_imp(quad_sig, x_offset, y_offset, z_offset, x_bias, y_bias, z_bias, window, num_windows)
    % Simulates  the effect of a quadrature pulse with z readout with imperfect conditions.
    % Zero relaxation, using hydrogens, readout along z.
    % 
    % Quad_sig should be an angle in radians, which is the angle off of the x
    % axis. Pass "45" for default pi/4.
    %
    % x, y and z bias are the bias fields in the system, in microgauss
    % x, y and z offsets are amount that the pulses are wrong in radians.
    %
    % Window is the length of time during acquisition in ms
    %
    % Num_windows is how many x-y pairs we'll acquire for. Total
    % acquisition time is num_windows * window * 4.
    %
    % Default values are: pi/4, 0, 0, 0, 0uG, 0uG, 0uG, 100ms, 8
    %
    % Usage: [out, t] = simulate_pulses(quad_sig, x_offset, y_offset, z_offset, x_bias, y_bias, z_bias, window, num_windows)
    %
    % out = 2D vector reads out x, y and z components.
    % 
    % t = time vector in ms.


    if nargin == 0 || quad_sig == 45
        quad_sig = pi/4;
    end

    if nargin < 2
        x_bias = 0;
    end

    if nargin < 3
        y_bias = 0;
    end

    if nargin < 4
        z_bias = 0;
    end

    if nargin < 5 || window < 1
        window = 100;
    end

    if nargin < 6 || num_windows < 1
        num_windows = 8;
    end
    
    if nargin < 7 || sr < 1
        sr = 500;
    end
  
    % Set up the x y and z coordinates.
    x = [1; 0; 0];
    y = [0; 1; 0];
    z = [0; 0; 1];
    
    % Set up the vector:
    v = z;
    
    % First pulse, pi/2 x to put it in the xy plane.
    v = RotVecArAxe(v, x, pi/2);
    
    % Second pulse, quad_sig pulse along z to give us the initial
    % quadrature signal.
    v = RotVecArAxe(v, z, quad_sig);
    
    % Now we put it into readout mode by tipping it to -z:
    v = RotVecArAxe(v, x, pi/2);
    
    % Set up the bias field
    bias_field = zeros(3, 1);
    bias_field(1) = x_bias;
    bias_field(2) = y_bias;
    bias_field(3) = z_bias;
    
    % Hydrogen gyromagnetic ratio is
    % 2.675e-5 rad * 1/(ms*uG)
    % We want our answer in radians.
    
    gam = 2.675e-5;
    mag = norm(bias_field); % Bias field magnitude in uG
    if mag == 0
        bias_field = z;
    else
        bias_field = bias_field./mag;
    end
    
    angle = mag*gam*window/sr; % Angle around bias_field traced out during window.
   
   % Initialize output and times vector.
    np = round(window*sr/1000); % How often to sample, round to nearest even integer.
    sr = np/window; % Adjust the sampling rate now that you've rounded.
    t = 1:(np*num_windows*4); % Time vector
    out = zeros(np*num_windows*4, 3); % Output vector.
    
    % Now loop through the quadrature signal, measuring as necessary,
    % skipping pulses.
    for i=0:4:num_windows*4
        for j = 1:np
            out(i*np+j, :) = v; % Magnetization
            v = RotVecArAxe(v, bias_field, angle); % Precession
        end
        
        v = RotVecArAxe(v, x, pi); % Apply the 180x pulse
        
        for j = 1:np
            out((i+1)*np+j, :) = v; % Magnetization
            v = RotVecArAxe(v, bias_field, angle); % Precession
        end

        v = RotVecArAxe(v, y, pi/2);    % Apply the 90y pulse
        
        for j = 1:np
            out((i+2)*np+j, :) = v; % Magnetization
            v = RotVecArAxe(v, bias_field, angle); % Precession
        end
        
        v = RotVecArAxe(v, x, pi); % Apply the last 180x pulse

        for j = 1:np
            out((i+3)*np+j, :) = v; % Magnetization
            v = RotVecArAxe(v, bias_field, angle); % Precession
        end
       
        % One more bias rotation, then the last pulse.
        v = RotVecArAxe(v, y, pi/2); % Apply the 90y pulse
    end


    
    
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
