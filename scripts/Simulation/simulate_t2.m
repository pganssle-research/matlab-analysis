function [out, t] = simulate_t2(t2, a_len, x_offset, y_offset, x_bias, y_bias, z_bias, np)
    % Simulates the t2 experiment.
    %
    % t2 = t2 of what you are looking at
    % a_len = acquisition length
    % np = number of points to acquire, default is 100.
    %
    % Offset is in us, bias is in uG, defaults are 0.
    %
    % We assume that the x and y pulse fields are 1.5G, z = 0.5G.
    %
    % Out is the t2 measurement output, t are the times.
    %
    % Usage: [out, t] = simulate_t2(t2, a_len, x_offset, y_offset, x_bias, y_bias, z_bias)

    if nargin == 0 || t2 <= 0
        t2 = 4.0;
    end
    
    if nargin < 2 || a_len <= 0
        a_len = 12.0;
    end
    
    if nargin < 3
        x_offset = 0;
    end

    if nargin < 4
        y_offset = 0;
    end

    if nargin < 5
        x_bias = 0;
    end

    if nargin < 6
        y_bias = 0;
    end

    if nargin < 7
        z_bias = 0;
    end
    
    if nargin < 8 || np <= 5
        np = 100;
    end
    
    % Sampling rate
    sr = np/a_len; % Sampling rate in Hz.

    % Set up the x y and z coordinates.
    x = [1; 0; 0];
    y = [0; 1; 0];
    z = [0; 0; 1];
    
    % Set up the offsets
    hgam = 2.6753e-2; % Gamma in MHz/G
    rx_off = 1.5*hgam*x_offset;
    ry_off = 1.5*hgam*y_offset;
    
    pi_x = pi+rx_off;
    pi2_x = pi/2 + rx_off/2;
    
    pi_y = pi+ry_off;
    pi2_y = pi/2 + ry_off/2;
    
    % Set up the bias field between pulses.
    bias_field = zeros(3, 1);
    bias_field(1) = x_bias;
    bias_field(2) = y_bias;
    bias_field(3) = z_bias;
    
    bf_zon = bias_field; %+ z*5e5; % Bias field during the pulse.
    
    mag = norm(bias_field); %Magnitude of the bias field by itself in microgauss.
    mag_zon = norm(bf_zon)/1e6; % Magnitude of the field in Gauss.
    
    % It doesn't like division by zero,
    % So normalize these on the condition that mag != 0.
    if mag == 0
        bias_field = z;
    else
        bias_field = bias_field/mag;
    end
    
    if mag_zon == 0
        bf_zon = z;
    else
        bf_zon = bf_zon/mag_zon;
    end
    
    % Time between pulses will be in increments of 100ms
    % Set up the angle it'll trace out between that time.
    z_ang = mag_zon*hgam*1000*1/(2*(sr*1000)); % hgam in kHz/G * sampling time = 1/2sr * 0.5G.
    b_ang = mag*hgam*1e-3; % Angle/time in kHz.
    
    % Loop through and generate the output vector.

    t = (0:(np-1))/sr; % Time vector in seconds.
    out = zeros(np, 1);
    
    for i=1:np
        v = z;
        v = RotVecArAxe(v, x, pi2_x); % Pi/2 x.
        v = RotVecArAxe(v, bf_zon, z_ang*i); % Evolution under the field.
        v = RotVecArAxe(v, x, pi_x); % Pi x.
        v = RotVecArAxe(v, y, pi_y); % Pi y.
        v = RotVecArAxe(v, x, pi_x); % Pi x.
        v = RotVecArAxe(v, bf_zon, z_ang*i); % Evolution under the field again.
        v = RotVecArAxe(v, x, pi2_x); % Pi/2 x.
        
        % Get the magnitude the same way we do.
        m1 = v(3); %-z
        v = RotVecArAxe(v, bias_field, b_ang*100); % Some bias field evolution between shots, 100ms.
        v = RotVecArAxe(v, x, pi_x); % Pi to get the +z;
        m2 = v(3); %+z
        
        m = m2-m1; % Magnetization at this point.
        m = m*exp(-(t(i)/t2)); % Decay.
        
        out(i) = m; % This is our output.
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
