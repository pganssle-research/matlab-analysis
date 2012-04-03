function [spectrum, f1, f2, fid, t1, t2] = sim_cosy(j_couplings, spins, enriched, pulse_length, bias_field, r)
% Simulates the 2D Zero Field "COSY" spectrum of a molecule.
%
% The inputs are as follows (n = number of spins):
% j_coupling = n x n matrix, the J-couplings in Hz
% 
% spins = nx1 vector telling me what nucleus each spin is.
%   0 = hydrogen
%   1 = 13-carbon
%   2 = 15-nitrogen
%   3 = 31-phosphorous
%   4 = 19-fluorine
%
% enriched is an nx1 vector telling me whether or not to use natural
%   abundance values for each nucleus. 0 = yes, 1 = 100% Defaults to 0 for
%   all
%
% pulse_length = how long the pulse should be, in units of hydrogen pi
%   pulses. so, for example a pi/2 hydrogen pulse should be passed 0.5.
%   If omitted, this defaults to 4.
%
% bias_field = 3x1 vector containing the bias fields in Gauss. Default is 0.
%
% r is the relaxation constant in Hz. Default is Inf.
%
% Usage:
% [spectrum, f1, f2, fid, t1, t2] = sim_cosy(j_couplings, spins, enriched, pulse_length, bias_field, r);

% Here are some gyromagnetic ratios:
hgam = 4257.6*2*pi; % Hydrogen gyromagnetic ratio in rad Hz/Gauss
cgam = 1070.5*2*pi; % 13-Carbon gyromagnetic ratio in rad Hz/Gauss
ngam = -431.6*2*pi; % 15-Nitrogen gyromagnetic ratio in rad Hz/Gauss
pgam = 1723.5*2*pi; % 31-Phosphorous gyromagnetic ratio in rad Hz/Gauss
fgam = 4005.3*2*pi; % 19-Fluorine gyromagnetic ratio in rad Hz/Gauss

% These are the natural abundances of 13-C and 15-N
c_ab = 0.01122; % 13-C abundance
n_ab = 0.003673; % 15-N natural abundance

% Establish parameters
h = 4.1367e-15; % h in eV * s
beta = 35.706; % 1/kb at 325K in eV^-1.
scale = h*beta; % hbar/beta at 325K in s

pi_pulse = 7.83e-5; % Hydrogen pi pulse length in seconds
pulse_strength = 1.5; % Pulse field strength in Gauss

% Check that the user gave us good data
s = size(j_couplings);
if numel(s) > 2
    fprintf('Incorrect size of j-coupling network matrix.\n');
    return
elseif s(2) ~= s(1)
    fprintf('J-coupling matrix must be square.\n');
    return
end

ns = s(1); % This is the number of spins

s = size(spins);

if s(2) > 1 && s(1) > 1
    fprintf('Incorrect shape of spins vector. It must be n x 1.\n');
    return
elseif s(2) > 1
    spins = spins';
elseif s(1) ~= ns
    fprintf('J coupling network has size %d. Spin vector has size %d. These must be the same.\n', ns, s(1));
    return
end

if nargin < 3
    enriched = zeros(ns, 1);
end

s = size(enriched);

if s(2) > 1 && s(1) > 1
    fprintf('Incorrect shape of abundance vector. It must be n x 1.');
    return
elseif s(2) > 1
    enriched = enriched';
    s = s';
end

if s(1) > ns
    fprintf('J coupling network and spin vector have size %d. Enriched vector has size %d. This should not exceed the number of spins.\n', ns, s(1));
    return
elseif s(1) < ns
    if s(1) == 0
        s(1) = 1; % Don't freak out if it's empty.
    end
    enriched(s(1):ns) = 0; % If it's not there, it's not enriched
end


if nargin < 4 || pulse_length == 0
    pulse_length = 4;
end

if nargin < 5 || numel(bias_field) ~= 3
    bias_field = [0 0 0];
end

if nargin < 6
    r = Inf;
end

pulse_length = pulse_length*pi_pulse; % Put this in the right units

% J coupling network, in Hz

% This one's for ethanol
if(0)

end

% Convert the j-coupling network into the right units
j_couplings = j_couplings*pi;

% Set up the abundances
abundances = zeros(ns, 1);
for n=1:ns
    if spins(n) == 1 && enriched(n) == 1
        abundances(n) = c_ab;
    elseif spins(n) == 2 && enriched(n) == 1
        abundances(n) = n_ab;
    else
        abundances(n) = 1;
    end
end


% Pre-allocate the higher-dimensional matrices
x_matrices = zeros(2^ns, 2^ns, ns);
y_matrices = zeros(2^ns, 2^ns, ns);
z_matrices = zeros(2^ns, 2^ns, ns);

% Set up the Pauli matrices
sig_x = (1/2)*[0 1; 1 0];
sig_y = (1/2)*1i*[0 -1; 1 0];
sig_z = (1/2)*[1 0; 0 -1];
id = eye(2);

% This generates the matrixes (xyz)_matrices
% z_matrices(:, :, n) is the sig_z matrix for the nth spin
for n = 1:ns
    if n == 1
        current_x = sig_x;
        current_y = sig_y;
        current_z = sig_z;
    else
        current_x = id;
        current_y = id;
        current_z = id;
    end
    
    if ns > 1        
        for m = 2:ns
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
    end
    
    x_matrices(:, :, n) = current_x;
    y_matrices(:, :, n) = current_y;
    z_matrices(:, :, n) = current_z;
end

% Sometimes you want to manipulate all of the spins in bulk.
% Generate these operators.
ip = unique(spins); % Listing of the isotopes present
ni = numel(ip); % Number of isotopes present

all_x = zeros(2^ns, 2^ns, ni);
all_y = zeros(2^ns, 2^ns, ni);
all_z = zeros(2^ns, 2^ns, ni);

for n = 1:ns
    for m = 1:ni
        if spins(n) == ip(m)
            all_x(:, :, m) = all_x(:, :, m) + x_matrices(:, :, n);
            all_y(:, :, m) = all_y(:, :, m) + y_matrices(:, :, n);
            all_z(:, :, m) = all_z(:, :, m) + z_matrices(:, :, n);
        end
    end
end
    
% Initial density matrix, assuming we start with magnetization along z
h_mag = (exp(-(scale*hgam*1e4)) - 1)/(exp(-scale*hgam*1e4)+1); % Magnetization of the hydrogens in the 1T field
c_mag = (exp(-(scale*cgam*1e4))-1)/(exp(-scale*cgam*1e4)+1); % Magnetization of the carbons in 1T field
n_mag = (exp(-(scale*ngam*1e4))-1)/(exp(-scale*ngam*1e4)+1); % Magnetization of the carbons in 1T field
p_mag = (exp(-(scale*pgam*1e4))-1)/(exp(-scale*pgam*1e4)+1); % Magnetization of the carbons in 1T field
f_mag = (exp(-(scale*fgam*1e4))-1)/(exp(-scale*fgam*1e4)+1); % Magnetization of the carbons in 1T field


h_frac = 1; % Set hydrogen magnetization to 1.
c_frac = c_mag/h_mag; % Relative proportion of carbon magnetization
n_frac = n_mag/h_mag; % Relative proportion of nitrogen magnetization
p_frac = p_mag/h_mag; % Relative proprtion of phosphorous magnetization
f_frac = f_mag/h_mag; % Relative proportion of fluorine magnetization


% First generate the hamiltonians

% Generate the J coupling portion of the Hamiltonian
h_j = zeros(2^ns); % Pre-allocate the hamiltonian
for n = 1:ns
    for m = 1:ns
        h_j = h_j + j_couplings(n, m)*dot_product(x_matrices, y_matrices, z_matrices, n, m); % There may be some better way to do this than nested for loops
    end
end

% Generate the pulse and free evolution Hamiltonians
h_pulse = zeros(2^ns); % Hamiltonian during the pulse
h_off = zeros(2^ns); % Hamiltonian when the pulse is off
readout_matrix = zeros(2^ns); % This is for readout later, not a Hamiltonian
for m = 1:numel(ip)
    n = ip(m);
    if n == 0
        gam = hgam;
    elseif n == 1
        gam = cgam;
    elseif n == 2
        gam = ngam;
    elseif n == 3
        gam = pgam;
    elseif n == 4
        gam = fgam;
    end
    
    n = n+1; % Because matlab uses a 1-based index.
    
    h_off = h_off + gam*(all_x(:, :, n)*bias_field(1) + all_y(:, :,n)*bias_field(2) + all_z(:, :, n)*bias_field(3)); % Add in the effect of the bias fields
    h_pulse = h_pulse + all_x(:, :, n)*gam*pulse_strength; % Pulse is along x, this is arbitrary
    readout_matrix = readout_matrix + all_z(:, :, n)*gam; % Readout matrix
end

h_off = h_off + h_j; % Finally the all-important J-coupling portion
h_pulse = h_pulse + h_off; % All the stuff that's present during h_off is also present during the pulse. This should be a perturbation.

% Get some information about the acquisition we'll need
np1 = 256; % Number of points in the direct dimension
np2 = 256; % Number of points in the indirect dimension

sr1 = 700; % Sampling rate in the direct dimension in Hz
sr2 = 700; % Sampling rate in the indirect dimension in Hz

t1 = (1:np1)/sr1; % Direct time vector
t2 = (1:np2)/sr2; % Indirect time vector

t1_int = 1/sr1; % Interval of time between direct measurements in seconds
t2_int = 1/sr2; % Interval of time between indirect measurements in seconds

% Now generate the propagators
pulse = expm(-1i*h_pulse*pulse_length); % Time = 4 pi on hydrogen
zf1 = expm(-1i*h_off*t1_int); % Time = 1/sr1. Operate on this between each direct dimension point.
zf2 = expm(-1i*h_off*t2_int); % Time = 1/sr2. Operate on this between each indirect dimension point

% Now we generate the FID

% Start out with the initial density matrix
% We want it to be a sum over the I_{iz} operators, weighted by the
% abundance and magnetized fraction.
rho = zeros(2^ns);
for n=1:ns
  if spins(n) == 0
    rho = rho + h_frac*abundances(n)*z_matrices(:, :, n);
  elseif spins(n) == 1
    rho = rho + c_frac*abundances(n)*z_matrices(:, :, n);
  elseif spins(n) == 2
    rho = rho + n_frac*abundances(n)*z_matrices(:, :, n);
  elseif spins(n) == 3
    rho = rho + p_frac*abundances(n)*z_matrices(:, :, n);
  elseif spins(n) == 4
    rho = rho + f_frac*abundances(n)*z_matrices(:, :, n);
  end
end

% Start off by applying the pulse to rho
rho2 = pulse*rho*pulse';

% This gives us the full FID by evolving the density matrix a bit at a
% time, then getting us the readout.
fid = zeros(np1, np2); % Pre-allocate the FID matrix
for m = 1:np2
    rho2 = zf2*rho2*zf2'; % Evolve for interval t in the indirect dimension
    
    rho1 = pulse*rho2*pulse'; % Starting from here, give us our second pulse
    
    for n = 1:np1
       rho1 = zf1*rho1*zf1';
       fid(n, m) = trace(readout_matrix*rho1);
    end
    
    % This is just to give us a little progress bar.
    if rem(m, ceil(np2/20)) == 0
        progress = (m/np2)*10;
        progress = floor(progress);
        fprintf('%d', progress);
    end
end
fprintf('\n');

relaxation = exp(-t1'/r)*exp(-t2/r);

fid = fid.*relaxation;

spec = fft2(fid); % Get the spectrum
f1 = 0:(sr1/(np1-1)):sr1; % Calculate the frequencies
f2 = 0:(sr2/(np2-1)):sr2;

% For convenience, I'm going to throw away the second half of the spectrum
% now. In the future, make sure these numbers of points are always powers
% of two so that this won't ever cause issues.

spectrum = spec(1:(np1/2), 1:(np2/2));
f1 = f1(1:(np1/2));
f2 = f2(1:(np2/2));

function o = dot_product(x, y, z, n, m)
% Generates the dot product between the two things you feed this.
x_n = squeeze(x(:, :, n));
x_m = squeeze(x(:, :, m));

y_n = squeeze(y(:, :, n));
y_m = squeeze(y(:, :, m));

z_n = squeeze(z(:, :, n));
z_m = squeeze(z(:, :, m));

o = x_n*x_m + y_n*y_m + z_n*z_m;