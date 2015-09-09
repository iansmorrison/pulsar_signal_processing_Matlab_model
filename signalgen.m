function signalgen()
% Generates a file containing dual noise vectors with phase-dependent
% partial polarization. File is 32-bit floating point with polarizations
% interleaved at each time step. 
 
% DATA SETTINGS
%
% fname     - Ouput filename
% hdrsize   - Header size
% hdrtype   - Data type for header ('uint8' = byte)
% ntype     - Data type for each element in a pair ('single' = float)
% Nout      - Length of each output vector
% nbins     - Number of bins within a pulse period
% npol      - Number of polarizations (should always be 2 when calc Stokes)
% nseries   - Number of forward FFT's to perform
% noise     - Set to 0.0 for no noise, 1.0 for noise (max(S/N)=1)
% dformat   - Specifies conversion TO real or complex data
% shift     - Select whether an fftshift is used before the inverse FFT
%             (don't shift if PFB is in the signal chain)
%
% INSTRUMENT SETTINGS
% f0        - Centre frequency (MHz) 
% f_sample_out - Sampling frequency of output data (MHz)
% 
% PULSAR SETTINGS
% Dconst    - Dispersion constant, s.MHz^2/(pc/cm^3)
% DM        - Dispersion measure, pc/cm^3
% pcal      - Pulsar period (s) and other params in a structure
% t0        - Absolute time (in seconds) of first time element 
%
% OUTPUTS:
% --------
%
%    fname -  file containing two interleaved floating point test vectors
%
% Description:
% ------------
% Generates a file containing dual noise vectors with phase-dependent
% partial polarization. File is 32-bit floating point with polarizations
% interleaved at each time step. 
% 
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% D. Hicks         04-Jul-2014  Original version
% I. Morrison      31-Jul-2015  Added noise parameter
%                               Added optional fftshift before inverse FFT
% ----------------------------------------------------------------------
 
%=============
 
fname = 'simulated_pulsar.dump';

hdrsize = 4096; % Header size
hdrtype = 'uint8'; % Data type for header ('uint8' = byte)
ntype = 'single'; % Data type for each element in a pair ('single' = float)
Nout = 2^20; %Length of each output vector
nbins = 2^10; % Number of bins within a pulse period
npol = 2; % Number of polarizations (should always be 2 when calc Stokes)
nseries = 30; % Number of FFT's to perform
noise = 0.0;  % 0.0 for no noise, 1.0 for noise (max(S/N)=1)
dformat = 'complextoreal'; %specifies conversion TO real or complex data
%dformat = 'complextocomplex'; %specifies conversion TO real or complex data
shift = 0; % performs an fftshift before the inverse FFT
 
% Instrument settings
f0 = 1405.; % Centre frequency (MHz)

% Set bandwidth to 8 x 10 MHz, for testing with 8-channel channelizer
f_sample_out = 80.; % Sampling frequency of output (MHz)

% Pulsar settings
Dconst = 4.148804E3; % s.MHz^2/(pc/cm^3)
DM = 2.64476; % pc/cm^3
pcal = struct('a',0.00575745,'b',0.0);% Pulsar period (s) and other params
t0 = 0.0; % Absolute time (in seconds) of first time element 
 
%===============
%Multiplying factor going from input to output type
switch dformat
    case 'complextoreal'
        Nmul = 2; 
    case 'complextocomplex'
        Nmul = 1;
    otherwise
        warning('Conversion should be complextoreal or complextocomplex.');
end;
 
Tout = 1/abs(f_sample_out)*1E-6; % Sample spacing of output (seconds)
df = f_sample_out/Nmul; % Bandwidth/Nyquist frequency (MHz)
Tin = Tout*Nmul; % Time spacing between input data elements
Nin = Nout/Nmul; % Number of data elements in input time series
Pmul = 1/Nmul; % Power multiplication factor for all but the DC channel
 
%===============
% Create the dispersion kernel and determine the number of elements to be
% clipped off the beginning and end.
frange = [-df/2, df/2] + f0
 
% Get matrix to perform dispersion on complex array
[H, ~, n_hi, n_lo] = ...
         dispnmatrix(frange, abs(df), Nin, 1, Dconst*DM, Tin, sign(df));
 
% Calculate the number of elements in the clipped input array
nclip_in = Nin - n_lo - n_hi; 
% Calculate number of elements in the clipped output array
nclip_out = Nout - n_lo*Nmul - n_hi*Nmul; 
 
frac_lost = (n_lo + n_hi)/Nin; % fraction of array that's lost
fprintf('Lost fraction of time series = %f\n', frac_lost);
fprintf('Time series length = %f s\n', nclip_in*Tin);

%===============
% Calculate phase-dependent Stokes parameters and coherency matrix
% using the rotating vector model
[~, J] = rotvecmod(nbins,noise);
 
% Vector of relative times
trel = (0:Nin-1)*Tin;
 
%===============
% Open file for writing
fid = fopen(fname, 'w');

% Write header
hdr = zeros(hdrsize,1);
fwrite(fid, hdr, hdrtype);
 
for ii = 1:nseries,
    % Print loop number
    fprintf('Loop # %i of %i\n', ii, nseries);
    
    % Time vector
    if ii == 1,
        tt = t0 - n_hi*Tin + trel;
    else
        tt = ttclip(end) - (n_hi-1)*Tin + trel;
    end;
    
    tindex = findphase(tt, nbins, pcal); 
    index = unique(tindex);
    
    % Initialize data vector for this series
    z = zeros(Nin, npol, 'single');
    %iL = 1; %Starting index when looping through phases
    
    % Loop through groups of data that share the same phase. Random data 
    % in each group are generated from the same coherency matrix
    
    for jj = 1:length(index),
        %Get coherency matrix for this pulsar phase
        Jcoh = [J(index(jj),1), J(index(jj),3); ...
                J(index(jj),2), J(index(jj),4)];
        
        % Indices of elements with a given phase
        iphase = find(tindex == index(jj));
        nL = length(iphase);
        
        %Generate two randomly-phased, unit-length phasors  
        %z0 = exp(complex(0,1)*2*pi()*rand(nL,npol));
        z0 = sqrt(0.5)*[complex(randn(nL,1),randn(nL,1)), ...
                        complex(randn(nL,1),randn(nL,1))];
 
        %Generate covariant vectors via Cholesky decomposition
        zjj = z0*chol(Jcoh, 'upper');
        %z = transpose(chol(Jcoh, 'lower')*transpose(z0)); %alternative
        
        % Concatenate with data from other phases
        z(iphase, :) = zjj;
        %iL = iL + nL; % increment to next starting index in z
    end;
    
    % Forward FFT
    f1a = fft(z(:,1), Nin);
    f2a = fft(z(:,2), Nin);
    
    % Element-wise multiplication by dispersion matrix.
    f1a = f1a .* H;
    f2a = f2a .* H;
    
    % If complextoreal, then create a Hermitian array
    switch dformat
        case 'complextoreal'
            %Create Hermitian vector
            f1 = [real(f1a(1)); f1a(2:Nin)*Pmul; ...
                  imag(f1a(1)); flipud(conj(f1a(2:Nin)))*Pmul];
            f2 = [real(f2a(1)); f2a(2:Nin)*Pmul; ...
                  imag(f2a(1)); flipud(conj(f2a(2:Nin)))*Pmul]; 
        otherwise
            f1 = f1a;
            f2 = f2a;
    end;
    
    % Inverse FFT
    % Optionally include an fftshift before the inverse FFT, as needed
    if shift == 1,
        f1 = fftshift(f1);
        f2 = fftshift(f2);
    end;
    z1 = ifft(f1, Nout);
    z2 = ifft(f2, Nout);
    
    % Remove convolution overlap region
    ttclip = tt(1+n_hi : Nin-n_lo);
    z1clip = z1(1+n_hi*Nmul : Nout-n_lo*Nmul);
    z2clip = z2(1+n_hi*Nmul : Nout-n_lo*Nmul);
 
    % Interleave polarizations into a single vector
    switch dformat
        case 'complextoreal'    
            z = [z1clip, z2clip];
            dat = reshape(transpose(z),npol*nclip_out,1);
        case 'complextocomplex'
            z = [real(z1clip), imag(z1clip), real(z2clip), imag(z2clip)];
            dat = reshape(transpose(z),2*npol*nclip_out,1);
    end
    
    %Write vector to file
    fwrite(fid, dat, ntype);
end;
 
fclose(fid);
return
 
end



function [S, J, p] = rotvecmod(N, noise, showplot)
% Rotating vector model for pulsar emission
 
if ~exist('N','var'),
    N = 1024;
end;
 
esig = 5. ; % emission half angle (polar angle, degrees)
epeak = 0. ; % emission peak angle (polar angle, degrees)
flin = 0.3; % linear polarized fraction amplitude
 
zeta = 30.; % observing angle (degrees) relative to rotation axis
alpha = 40.; % magnetic axis (degrees) relative to rotation axis
 
pmin = -180.;
pmax = 180.;
 
% Angle of rotation: p=0 for aligned dipole. 
% This is equivalent to pulsar longitude or phase
p = transpose(linspace(pmin, pmax, N)); 
 
% Polarization angle w.r.t projected rotation axis from observing direction
%psi = atand(sind(alpha)*sind(p)./(sind(zeta)*cosd(alpha) - ...
%    sind(alpha)*cosd(zeta)*cosd(p)));
psi = atan2d(sind(alpha)*sind(p),  ...
    (sind(zeta)*cosd(alpha) - sind(alpha)*cosd(zeta)*cosd(p)));
 
% Polar observation angle in magnetic axis reference frame
cosO = cosd(p)*sind(zeta)*sind(alpha) + cosd(alpha)*cosd(zeta);
tanO = sqrt(1./(cosO.^2)-1);
 
% Polar emission angle in magnetic axis reference frame
thetaE = atand(1.5*(sqrt(1+(8/9)*tanO.^2) - 1)./tanO);
%thetaE = atand(1.5*(-sqrt(1+(8/9)*tanO.^2) - 1)./tanO);
 
% Intensity (model-based assumption)
S0 = (1./sqrt(2*pi()*esig^2))*exp(-(thetaE-epeak).^2/(2.*esig^2));
S0 = S0/max(S0); %normalize max to 1
 
% Linear polarization fraction (model-based assumption)
L = flin*S0.*cosd(thetaE);
 
% Other Stokes parameters
S1 = L.*cosd(2*psi);
S2 = L.*sind(2*psi);
S3 = -(1-flin)*S1; % Fake circular polarization to avoid zero signal
%S3 = single(zeros(N,1)); % Zero circular polarization component
 
% Add noise, typically such that max(S/N) = 1
S0 = S0 + noise;
 
% Normalize Stokes 4-vector so that S0 = 1. 
factor = max(S0);
S0 = S0/factor;
S1 = S1/factor;
S2 = S2/factor;
S3 = S3/factor;
 
% Create Coherency matrix
Jxx = 0.5*(S0 + S1);
Jyy = 0.5*(S0 - S1);
Jxy = 0.5*(S2 + 1i*S3);
Jyx = 0.5*(S2 - 1i*S3);
 
% Plot results, if requested. Useful for debugging.
if exist('showplot','var'),
    clf();
 
    subplot(2,2,1);
    plot(p, transpose([S0, S1, S2, S3])); 
    legend('S0', 'S1', 'S2', 'S3');
    xlabel('Longitude (degrees)','FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
 
    subplot(2,2,3);
    plot(S0, transpose([S1, S2]));
    hleg1 = legend('S1', 'S2');
    set(hleg1,'Location','NorthWest')
    axis([0, 2, -Inf, Inf]);
    xlabel('S0','FontSize', 12, 'FontWeight', 'bold');
    ylabel('S1 or S2','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
 
    subplot(2,2,2);
    plot(p, transpose([Jxx, Jyy, real(Jxy), imag(Jxy)]));
    legend('Jxx', 'Jyy', 'Real(Jxy)', 'Imag(Jxy)');
    xlabel('Longitude (degrees)','FontSize', 12, 'FontWeight', 'bold');
    ylabel('Amplitude','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
 
    subplot(2,2,4);
    plot(S1, S2, 'b');
    xlabel('S1','FontSize', 12, 'FontWeight', 'bold');
    ylabel('S2','FontSize', 12, 'FontWeight', 'bold');
    set(gca,'FontSize', 12, 'FontWeight', 'bold');
end;
 
S = [S0, S1, S2, S3];
J = [Jxx, Jyx, Jxy, Jyy];
return
end
