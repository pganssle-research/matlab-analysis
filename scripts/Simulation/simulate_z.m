function [out_x out_y, t] = simulate_z(np, sr, field, bias_vec, num_wins, window, start_win)
    % Simulates an indirect acquisition with 1D quadrature readout.
    % 180s are along x.
    %
    % np = number of points, default is 20
    % sr = sampling rate in kHz, default is 20kHz
    % field = z field in Gauss, default is 0.5G
    % bias_vec = field bias, [x;y;z], in uG
    % num_wins = the number of windows to acquire, default is 1
    % start_win = the window to start with, default is 1
    % window = how wide is the window, in ms, default is 100
    %
    % t = time vector
    % out_x = (np x num_wins) vector of magnetization
    % out_y = (np x num_wins) vector of magnetization
    %
    % Usage: [out_x out_y t] = simulate_z(np, sr, field, bias_vec, num_wins, start_win, window);
    
    % Set up the defaults:
    if nargin == 0
        np = 20;
    end
    
    if nargin < 2
        sr = 20;
    end
    
    if nargin < 3
        field = 0.5;
    end
    
    if nargin < 4
        bias_vec = [0; 0; 0];
    end
    
    if nargin < 5
        num_wins = 1;
    end
    
    if nargin < 7
        start_win = 1;
    end
    
    if nargin < 6
        window = 100;
    end
            
  
    % Calculate the z angle for hydrogens.
    gam = 26.75; % Hydrogen gyromagnetic ratio in rad * kHz/G.
    len = 1/sr; % Length of time between acquisitions in kHz.
    
    angle = len*gam*field; % Angle in radians.
    
    % Initialize outputs:
    out_x = zeros(np, num_wins);
    out_y = zeros(np, num_wins);
    t = (1:np)/sr; % Time vector
    
    for i = 1:np
        % Out will be a 2D vector where the rows are the windows. 
        % out(:, 1) = x component.
        % out(:, 2) = y component.    
        out = simulate_pulses(angle*i, bias_vec(1), bias_vec(2), bias_vec(3), window, start_win+num_wins);
        
        out = out(start_win:start_win+num_wins-1, :);
                
        out_x(i, :) = out(:, 1);
        out_y(i, :) = out(:, 2);
    end    
    