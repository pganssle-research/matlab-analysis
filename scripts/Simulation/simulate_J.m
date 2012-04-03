function [spectrum, f, fid, t] = simulate_J(j_couplings, spins, enriched, pulse_length, bias_field, r)
% Simulates the J spectrum of a molecule.
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
%   abundance values for each nucleus. 0 = This nucleus is not included in
%   the spectrum. 1 = This is enriched and should be included 100%
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
% [spectrum, f, fid, t] = simulate_J(j_couplings, spins, enriched, pulse_length, bias_field, r);

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
    s = fliplr(s);
end

if s(1) > ns
    fprintf('J coupling network and spin vector have size %d. Enriched vector has size %d. This should not exceed the number of spins.\n', ns, s(1));
    return
elseif s(1) < ns
    if s(1) == 0
        s(1) = 1; % Don't freak out if it's empty.
    end
    enriched(s(1)+1:ns) = 0; % If it's not there, it's not enriched
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

% Convert the j-coupling network into the right units
j_couplings = j_couplings*pi;

% Set up the abundances
% This is tricky, because if they pass "2" as any of the arguments, we need
% to run the simulation at least twice, then take the weighted linear combination of
% the two.
% If they pass '0' to a carbon or nitrogen, we need to pull those indices
% out of everything that uses them.
abundances = zeros(ns, 1);
k = 1;

for n=1:ns
    if enriched(n) == 0 && (spins(n) == 1 || spins(n) == 2) 
        R(k) = n;
        k = k+1;
    else
        abundances(n) = 1;
    end
end

k = k-1;

% If there are spins we need to ignore, pull them out
if exist('R', 'var')
    ns = ns-k;
    for m = 1:k
        abundances(R(m)) = [];
        j_couplings(R(m), :) = [];
        j_couplings(:, R(m)) = [];
        spins(R(m)) = [];
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

[h_jbs h_j2] = eig(h_j);

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
h_pulse = h_pulse;% + h_off; % All the stuff that's present during h_off is also present during the pulse. This should be a perturbation.

% Get some information about the acquisition we'll need
np = 4096; % Number of points in the direct dimension
sr = 1400; % Sampling rate in the direct dimension in Hz
t = (1:np)/sr; % Time vector
t_interval = 1/sr; % Interval of time between measurements in seconds
fid = zeros(np, 1); % Pre-allocate the FID matrix

% Now generate the propagators
pulse = expm(-1i*h_pulse*pulse_length); % Time = 4 pi on hydrogen
zf = expm(-1i*h_off*t_interval); % Time = 1/sr. Operate on this between each point.

% Now we generate the FID

% Start out with the initial density matrix
rho = zeros(2^ns); % Pre-allocate rho.

% We want it to be a sum over the I_{iz} operators, weighted by the
% abundance and magnetized fraction.
for n=1:ns
  if spins(n) == 0
    rho = rho + z_matrices(:, :, n);
  elseif spins(n) == 1
    rho = rho + z_matrices(:, :, n);
  elseif spins(n) == 2
    rho = rho + z_matrices(:, :, n);
  elseif spins(n) == 3
    rho = rho + z_matrices(:, :, n);
  elseif spins(n) == 4
    rho = rho + z_matrices(:, :, n);
  end
end


% Start off by applying the pulse to rho

rho = pulse*rho*pulse';
% rho = h_jbs'*rho*h_jbs;

if(0)
rho(1:3, 4:end) = 0;
rho(4:end, 1:3) = 0;

rho(1, 2:3) = 0;
rho(2:3, 1) = 0;
rho(2, 3) = 0;
rho(3, 2) = 0;

rho(12, 13:end) = 0;
rho(13:end, 12) = 0;
rho(13, 14:end) = 0;
rho(14:end, 13) = 0;
rho(14, 15:end) = 0;
rho(15:end, 14) = 0;
rho(15, 16) = 0;
rho(16, 15) = 0;

rho(11:16, 4:11) = 0;
rho(4:11, 11:16) = 0;
end

m1 = zeros(16) + eye(16);
%m1(1:3, 12:16) = 1;
%m1(12:16, 1:3) = 1;

%%m1(1:2, 16) = 1; % Mostly 1
%%m1(16, 1:2) = 1;

%%m1(1:2, 12) = 1; % Mostly 2
%%m1(12, 1:2) = 1;

%%m1(3, 14) = 1;
%%m1(14, 3) = 1;

%m1(3, 12) = 1;
%m1(12, 3) = 1;

%m1(1:3, [12 14 16]) = 1;
%m1([12 14 16], 1:3) = 1;

%m1(1:3, 1:3) = 1;
%m1(12:16, 12:16) = 1;
%m1(4:5, 4:5) = 1;

%m1(6:11, 6:11) = 1;

%m1(4:5, 12:16) = 1;
%m1(12:16, 4:5) = 1;

%%m1(1:3, 4:5) = 1;
%%m1(4:5, 1:3) = 1;

%%m1(6:11, 12:16) = 1;
%%m1(12:16, 6:11) = 1;

m1(5, 7:10) = 1; % Mostly 7
m1(7:10, 5) = 1;
m1(4, [7 8 11]) = 1; % Mostly 11
m1([7 8 11], 4) = 1;

%%m1(1:3, 6:11) = 1;
%%m1(6:11, 1:3) = 1;
 
%rho = rho.*m1;

%rho = h_jbs*rho*h_jbs';


% Get the readout matrix

% This gives us the full FID by evolving the density matrix a bit at a
% time, then getting us the readout.
for n = 1:np
    rho = zf*rho*zf'; % Evolve for interval t
    
    fid(n) = trace(readout_matrix*rho); % This is the FID trace
    
    % This is just to give us a little progress bar.
    if rem(n, ceil(np/25)) == 0
        if (n/np) < 0.25
            fprintf('1');
        elseif (n/np) < 0.5
            fprintf('2');
        elseif (n/np) < 0.75
            fprintf('3');
        else
            fprintf('4');
        end
    end
end
fprintf('\n');


fid = fid.*exp(-t'/r);

spec = fft(fid); % Get the spectrum
f = 0:(sr/(np-1)):sr; % Calculate the frequencies

% For convenience, I'm going to throw away the second half of the spectrum
% now. In the future, make sure these numbers of points are always powers
% of two so that this won't ever cause issues.

spectrum = spec(1:(np/2));
f = f(1:(np/2));

function o = dot_product(x, y, z, n, m)
% Generates the dot product between the two things you feed this.
x_n = squeeze(x(:, :, n));
x_m = squeeze(x(:, :, m));

y_n = squeeze(y(:, :, n));
y_m = squeeze(y(:, :, m));

z_n = squeeze(z(:, :, n));
z_m = squeeze(z(:, :, m));

o = x_n*x_m + y_n*y_m + z_n*z_m;