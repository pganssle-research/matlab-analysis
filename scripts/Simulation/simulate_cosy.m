function [spectrum, f1, f2, fid, t1, t2] = simulate_cosy()
% Simulates a 2D zero field "COSY", which starts with magnetization along
% z, followed by a 4pi hydrogen pulse, followed by a delay with evolution
% at zero field, followed by another 4pi pulse.
%
% Currently configured for ethanol
%
% Outputs are a 2D spectrum and frequencies
% f1 = frequency in the direct dimension
% f2 = frequency in the indirect dimension
%
% Optional outputs are fid and time, which are 2D vectors as well
% t1 = time in the direct dimension
% t2 = time in the indirect dimension
%
% Usage:
% [spectrum, f1, f2, fid, t1, t2] = simulate_cosy();

% Establish parameters
hgam = 4257.6; % Hydrogen gyromagnetic ratio in Hz/Gauss
cgam = 1070.5; % Carbon gyromagnetic ratio in Hz/Gauss
h = 4.1367e-15; % h in eV * s
beta = 35.706; % 1/kb at 325K in eV^-1.
scale = h*beta; % hbar/beta at 325K in s

pi_pulse = 7.83e-5; % Pi pulse length in seconds
pulse_strength = 1.5; % Pulse field strength in Gauss
bias_field = [0; 0; 0]; % Bias field offset in Gauss
r = 3; % Relaxation time constant in seconds.

% Seven spins, establish the base operators
% The order is as follows:
% 1 = alpha carbon
% 2 = beta carbon
% 3, 4 = alpha hydrogens
% 5, 6, 7 = beta hydrogens

% Pre-allocate the higher-dimensional matrices
x_matrices = zeros(2^7, 2^7, 7);
y_matrices = zeros(2^7, 2^7, 7);
z_matrices = zeros(2^7, 2^7, 7);

% Set up the Pauli matrices
sig_x = (1/2)*[0 1; 1 0];
sig_y = (1/2)*i*[0 -1; 1 0];
sig_z = (1/2)*[1 0; 0 -1];
id = eye(2);

% This generates the matrixes (xyz)_matrices
% z_matrices(:, :, n) is the sig_z matrix for the nth spin
for n = 1:7
    if n == 1
        current_x = sig_x;
        current_y = sig_y;
        current_z = sig_z;
    else
        current_x = id;
        current_y = id;
        current_z = id;
    end
    
    for m = 2:7
        if m == n
            current_x = kron(current_x, sig_x);
            current_y = kron(current_y, sig_y);
            current_z = kron(current_z, sig_z);
        else
            current_x = kron(current_x, id);
            current_y = kron(current_y, id);
            current_z = kron(current_z, id);
        end
    end
    
    x_matrices(:, :, n) = current_x;
    y_matrices(:, :, n) = current_y;
    z_matrices(:, :, n) = current_z;
end

% Sometimes you want to manipulate all of the spins in bulk.
% Generate these operators.
all_carbons_x = squeeze(sum(x_matrices(:, :, 1:2), 3));
all_carbons_y = squeeze(sum(y_matrices(:, :, 1:2), 3));
all_carbons_z = squeeze(sum(z_matrices(:, :, 1:2), 3));

all_hydrogens_x = squeeze(sum(x_matrices(:, :, 3:end), 3));
all_hydrogens_y = squeeze(sum(y_matrices(:, :, 3:end), 3));
all_hydrogens_z = squeeze(sum(z_matrices(:, :, 3:end), 3));
    
% J coupling network, in Hz
j_couplings = [0 0 140.4 140.4 -4.6 -4.6 -4.6;
 0 0 -2.4 -2.4 125.2 125.2 125.2;
 140.4 -2.4 0 0 7.1 7.1 7.1;
 140.4 -2.4 0 0 7.1 7.1 7.1;
 -4.6 125.2 7.1 7.1 0 0 0;
 -4.6 125.2 7.1 7.1 0 0 0;
 -4.6 125.2 7.1 7.1 0 0 0];

% Initial density matrix, assuming we start with magnetization along z
c_mag = (exp(-(scale*cgam*1e4))-1)/(exp(-scale*cgam*1e4)+1); % Magnetization of the carbons in 1T field
h_mag = (exp(-(scale*hgam*1e4)) - 1)/(exp(-scale*hgam*1e4)+1); % Magnetization of the hydrogens in the 1T field

h_frac = 1; % Set hydrogen magnetization to 1.
c_frac = c_mag/h_mag; % Relative proportion of carbon magnetization

a_prop = 1; % 1 if this is alpha-labeled, 0.01109 for natural abundance
b_prop = 0.01109; % 1 if this is beta-labeled, 0.01109 for natural abundance

rho_init = h_frac*all_hydrogens_z + c_frac*a_prop*z_matrices(:, :, 1) + c_frac*b_prop*z_matrices(:, :, 2); % Initial rho is along z

% Generate the J coupling portion of the Hamiltonian
h_j = zeros(2^7, 2^7);
for n = 1:7
    for m = 1:7
        h_j = h_j + j_couplings(n, m)*dot_product(x_matrices, y_matrices, z_matrices, n, m);
    end
end

% Generate the pulse and free evolution Hamiltonians
h_pulse = all_carbons_y*cgam*pulse_strength + all_hydrogens_y*hgam*pulse_strength + h_j; % J coupling still present
h_off = (all_carbons_x*bias_field(1) + all_carbons_y*bias_field(2) + all_carbons_z*bias_field(3))*cgam;
h_off = h_off + (all_hydrogens_x*bias_field(1) + all_hydrogens_y*bias_field(2) + all_hydrogens_z*bias_field(3))*hgam;
h_off = h_off + h_j; % The first term is all the bias fields, the second term is the j coupling.

% Now generate the evolution matrices (I imagine this will be the long part
pulse = exp(-i*h_pulse*4*pi_pulse); % Time = 4 pi on hydrogen
pulse_t = pulse';

zf = exp(-i*h_off); % The times vary, but do this part, then raise it to the t power later.

% Calculate the eigenvectors to diagonialize the zf matrix
[zf_ev zf_evals] = eig(zf); % There will be roundoff error here.


% Now we generate the FIDs
np_dir = 512; % Number of points in the direct dimension
np_indir = 512; % Number of points in the indirect dimension

sr_dir = 700; % Sampling rate in the direct dimension in Hz
sr_indir = 700; % Sampling rate in the indirect dimension in Hz

t_vec_1 = (1:np_dir)/sr_dir;
t_vec_2 = (1:np_indir)/sr_indir;

fid = zeros(np_dir, np_indir); % Pre-allocate the FID matrix

% Get the zf in the diagonal basis.
zf_diag = zf_ev'*zf*zf_ev;
zf_diag = zf_diag.*eye(size(zf_diag)); % It should be diagonal up to a rounding error. So we don't exponentiate the rounding errors, I'm eliminating them.

% This gives us the full FID
for m = 1:np_indir
    %First calculate the evolution during the indirect time.
    zf_indir = zf_diag.^t_vec_2(m); % Do the exponentiation bit now.
    zf_indir = zf_ev*zf_indir*zf_ev'; % Bring it back to the original basis.
    
    rho_t2 = pulse*rho_init*pulse_t; % This is the density matrix immediately after the pulse.
    rho_t2 = zf_indir*rho_t2*zf_indir'; % Evolution in zero field.
    rho_t2 = rho_t2*exp(-t_vec_2(m)/r); % Relaxation/
    rho_t2 = pulse*rho_t2*pulse_t; % Density matrix going into the direct dimension.
    
    for n = 1:np_dir
        zf_dir = zf_diag.^t_vec_1(n); % Exponentiation in the diagonal basis
        zf_dir = zf_ev*zf_dir*zf_ev'; % Evolution in zero field
        
        rho_t1 = zf_dir*rho_t2*zf_dir'; % Evolution now
        rho_t1 = rho_t1*exp(-t_vec_1(n)/r); % Relaxation again.       
       
        % Finally, we need to get the trace of the magnetization vector.
        readout_matrix = all_hydrogens_z*hgam + all_carbons_z*cgam;
        fid(n, m) = trace(readout_matrix*rho_t1); % All that work for a single point in 25,000 or so
    end
    if rem(m, 20) == 0
        fprintf('|');
    end
end
fprintf('\n');

% FIDs are done, so now initialize the time vectors
t1 = t_vec_1;
t2 = t_vec_2;

spec = fft2(fid); % Get the spectrum
f1 = 0:(sr_dir/(np_dir-1)):sr_dir; % Calculate the frequencies
f2 = 0:(sr_indir/(np_indir-1)):sr_indir;

% For convenience, I'm going to throw away the second half of the spectrum
% now. In the future, make sure these numbers of points are always powers
% of two so that this won't ever cause issues.

spectrum = spec(1:(np_dir/2), 1:(np_indir/2));
f1 = f1(1:(np_dir/2));
f2 = f2(1:(np_indir/2));


function o = dot_product(x, y, z, n, m)
% Generates the dot product between the two things you feed this.
x_n = squeeze(x(:, :, n));
x_m = squeeze(x(:, :, m));

y_n = squeeze(y(:, :, n));
y_m = squeeze(y(:, :, m));

z_n = squeeze(z(:, :, n));
z_m = squeeze(z(:, :, m));

o = x_n*x_m + y_n*y_m + z_n*z_m;