function dat = signalproc(Nbits,remove_rfi)
%
% Given a stream of baseband data with 2 interleaved polarizations this 
% calculates the 4 Stokes parameters folded over a pulsar period. Data 
% are coherently de-dispersed over one or more frequency channels.  
%
% Inputs:
% -------
% (For now these are written directly into the procedure but they should
% ultimately be provided as inputs based on a defined control interface).
%
%   DATA SETTINGS
%
%   fname       - Input filename
%   hdrsize     - Size of header (in bytes)
%   hdrtype     - Header type (usually bytes = 'uint8')
%   Nin         - Number of elements in a given data read (power of 2)
%   ntype       - Data type (e.g. float = 'single')
%   nbytes      - Number of bytes per data point
%   npol        - Number of polarisations (should always be 2)
%   nbins       - Number of bins within a pulsar period
%   nfreq       - Number of frequency channels to divide the spectrum into
%   nseries     - Number of times the read loop should be run
%
%   INSTRUMENT SETTINGS
% 
%   dformat     - Format of incoming data (real or complex)
%   shift       - Select whether an fftshift is used after the FFT
%                 (shift needed if PFB is in the signal chain)
%   f0          - Centre frequency 
%   f_sample_in - Input data sampling frequency
%   resample    - Select whether to resample to the usable bandwidth after
%                 de-dispersion (for oversampled input)
%
%   RFI PARAMETERS
%
%   Mrfi        - Number of samples in the ensemble
%   pfn         - Probability of a false positive
%
%   PULSAR SETTINGS 
%
%   Dconst      - Dispersion constant in us.MHz^2/(pc/cm^3)
%   DM          - Dispersion measure in pc/cm^3
%   pcal        - Structure containing pulsar information
%   t0          - Absolute time (seconds) of first time element
%
% Outputs:
% --------
%
%    dat -  structure containing folded Stokes parameters and 
%           analysis metadata.
%
% Description:
% ------------
% Calculates the 4 Stokes parameters as a function of pulsar period. This
% is a Matlab implementation of the essential analysis components of DSPSR. 
% 
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% D. Hicks         21-Apr-2014  Original version
% D. Hicks         04-Jul-2014  Reads complex data, reduces bits, RFI
% I. Morrison      31-Jul-2015  Various provisions for oversampled input
%                               Added optional fftshift before inverse FFT
% ----------------------------------------------------------------------
 
 
close all; clear all; clc;

%=============
% Input parameters for data, instrument, and pulsar. For now these are
% written directly into the procedure instead of being passed into the 
% function as an input or read from the data header. This can be updated
% when the control interfaces are defined.
 
% Data settings
%fname = 'simulated_pulsar.dump';
fname = 'os_channelized_pulsar.dump';

hdrsize = 4096; % Header size
hdrtype = 'uint8'; % Data type for header ('uint8' = byte)
Nin = 2^18; % Number of elements per initial FFT
ntype = 'single'; % Data type for each element in a pair ('single' = float)
nbytes = 4; % Number of bytes per data point for ntype
npol = 2; % Number of polarisations (should always be 2 when calc Stokes)
nbins = 2^10; % Number of bins within a pulse period
nfreq = 1; % Number of frequency channels to divide the spectrum into
nseries = 5; % Number of time series to read and fft.
 
% RFI parameters
Mrfi = 2^5; % Number of samples in ensemble
pfn = 0.001; % Probability of false positive
 
% Instrument settings
%dformat='realtocomplex'; % Specifies conversion OF real or complex data
dformat='complextocomplex'; % Specifies conversion OF real or complex data
shift = 1; % 1 to perform an fftshift after the FFT, 0 if not
f0 = 1405.; % MHz
os_factor = 8.0/7.0; % Set to 1.0 for no oversampling
f_sample_in = 10.0*os_factor; % Sampling frequency of input (MHz)
% select part of the total Nyquist BW to be de-dispersed and used (MHz)
usable_bandwidth = 10.0; %abs(f_sample_in);
resample = 0;   % 1 to resample to the usable BW after de-dispersion,
                % 0 to stay oversampled
 
% Pulsar settings
Dconst = 4.148804E3; % s.MHz^2/(pc/cm^3)
DM     = 2.64476; % pc/cm^3
pcal   = struct('a',0.00575745, 'b', 0.0); % Pulsar period and other params
t0     = 0.; % Absolute time (in seconds) of first time element 
 
%Time offset (in seconds) at which to start file read (for debugging)
%tau = 0.4; 
 
%=============
% Default value of Nbits is 32 (32-bit floating point)
if ~exist('Nbits','var'),
    Nbits = 32;
end;
%=============
% Default is to not remove RFI
if ~exist('remove_rfi','var'),
    remove_rfi = 0;
end;
 
% Calculate spectral Kurtosis cutoffs
if remove_rfi == 1,
    % find SK cutoffs
    [Klo, Khi] = RFI_getSKlims(Mrfi,pfn); 
    % fast approximation of cutoff 
    %Klo = 1-sqrt(4/M)*3; Khi = 1+sqrt(4/M)*3; 
end;        
 
%=============
% Set up parameters depending on whether incoming data is real or complex
switch dformat
    case 'realtocomplex' 
        Nmul = 2; % Multiplying factor in converting from real to complex
        NperNin = 1; % 1 data point per polariz per incoming time step
    case 'complextocomplex'
        Nmul = 1;
        NperNin = 2; % 1 real + 1 imag point per pol per incoming time step 
    otherwise
        warning('Conversion should be realtocomplex or complextocomplex.');
end
 
Tin  = 1/abs(f_sample_in)*1E-6; % Sample spacing of input (seconds)
df   = f_sample_in/Nmul; % Bandwidth/Nyquist frequency (MHz)
Tout = Tin*nfreq*Nmul; % Time spacing between output data elements
Nout = Nin/nfreq/Nmul; % Number of data elements in output time series
Pmul = Nmul; % Power multiplication factor for all but the DC channel

%=============
% Create de-dispersion kernel with requested number of frequency channels
% Determine number of elements to be clipped off beginning and end.
frange = [-df/2, df/2] + f0

% Get matrix to perform de-dispersion on complex array
[H, f, n_hi, n_lo] = dispnmatrix(...
        frange, usable_bandwidth, Nout, nfreq, Dconst*DM, Tout, -sign(df));
 
% Length of output after clipping                     
nclip = Nout - n_lo - n_hi;
    
% Parameters after re-sampling (if done)
if resample == 1,
    bw = f_sample_in;
    disc = (abs(bw) - usable_bandwidth)/abs(bw)/2.0;
    newNout = ceil(Nout*(1.0-2.0*disc));
    newTout = (1/abs(usable_bandwidth)*1E-6)*nfreq*Nmul;
    % retain same clipping in time
    nclip = round(nclip/(abs(f_sample_in)/usable_bandwidth));
    % may need to manually tweak new_n_lo & new_n_hi for required tt length
    new_n_lo = round(n_lo/os_factor + 0.5);
    new_n_hi = ceil(n_hi/os_factor + 0.5);
else
    newNout = Nout;
    newTout = Tout;
    new_n_lo = n_lo;
    new_n_hi = n_hi;
end;

frac_lost = (n_lo + n_hi)/Nout; % Fraction of array that's lost
fprintf('Lost fraction of time series = %f\n', frac_lost);
fprintf('Time series length = %f s\n', nclip*newTout);

%==============
fid = fopen(fname);
 
% Read header
fread(fid, hdrsize, hdrtype);
%disp(transpose(native2unicode(hdr))); % Show header
 
% Shift starting point for reading file (useful for debugging)
%pt = round(tau/Tin)*npol*nbytes*NperNin;
%fseek(fid, pt, 'cof');

% Sum of folded Stokes data 
Vsum = zeros(nbins, nfreq, 4, 'single'); 
% Sum of squares of folded Stokes data
Vsum2 = zeros(nbins, nfreq, 4, 'single');
% Number of points per folded bin
Vn   = zeros(nbins,1, 'single'); 
% Vector of times relative to first element
trel = (0:nclip-1)*newTout;

%==============
% Main Loop
for ii = 1:nseries,
    % Print loop number
    fprintf('Loop # %i of %i\n', ii, nseries);
    
    if ii == 1,
        % Vector of times for the first (clipped) array to be read
        tt = t0 + n_hi*Tout + trel; %seconds
    else
        % Vector of times for the next (clipped) array to be read
        tt = tt(end) + newTout + trel;
        % Shift file pointer back by # of clipped elements for next read
        fseek(fid, -(n_hi+n_lo)*(Nin/Nout)*npol*nbytes*NperNin, 'cof');
    end;
        
    % Read stream of voltages, convert to 32-bit for speed
    % This forms a single column
    Vstream = single(fread(fid, npol*Nin*NperNin, ntype));
    
    if feof(fid)
        error('Error - hit end of input file!');
    end;
    
    %====================
    %Degrade to a given number of bits
    
    if Nbits < 32,
        Vstream = reducebits(Vstream,Nbits);
    end;
    %====================
    % Parse real and imag components if incoming data is complex
    
    switch dformat
        case 'complextocomplex'
            Vstream = reshape(Vstream, 2, []);
            Vstream = complex(Vstream(1,:), Vstream(2,:));
    end;
    
    % Separate data into different polarisations: Vdat(1,:) and Vdat(2,:)
    Vdat = reshape(Vstream, npol, []);
    
    %====================
    % Remove RFI
    if remove_rfi == 1,
        Nrfi = Nin/Mrfi; % Number of elements in a given sample
        % Remove channels that are outside the SK cutoffs. 
        [z1, Nlost1, SK1] = RFI_remove(Vdat(1,:),Nrfi,Mrfi,Klo,Khi);
        [z2, Nlost2, SK2] = RFI_remove(Vdat(2,:),Nrfi,Mrfi,Klo,Khi);
        Vdat = [z1; z2];
        fprintf('Mean SK1 = %f\n', mean(SK1));
        fprintf('Mean SK2 = %f\n', mean(SK2));
        fprintf('Fraction of lost channels in pol 1 = %f\n', Nlost1/Nrfi);
        fprintf('Fraction of lost channels in pol 2 = %f\n', Nlost2/Nrfi);
    end;
    %====================
    % Complex FFT for each polarisation; break & downshift into nfreqs.
    % Coherently de-disperse then inverse fft to get the analytic signal
    
    % Form analytic signal by taking the first half and multiplying by Pmul
    switch dformat
        case 'realtocomplex'
            F1all = fft(Vdat(1,:), Nin);
            F2all = fft(Vdat(2,:), Nin);
            F1 = [complex(F1all(1), F1all(Nout*nfreq+1)), ...
                F1all(2:Nout*nfreq)*Pmul];
            F2 = [complex(F2all(1), F2all(Nout*nfreq+1)), ...
                F2all(2:Nout*nfreq)*Pmul];
        case 'complextocomplex'
            F1 = fft(Vdat(1,:), Nin);
            F2 = fft(Vdat(2,:), Nin);
            % Optionally fftshift after the FFT, as needed
            if shift == 1,
                F1 = fftshift(F1);
                F2 = fftshift(F2);
            end;
    end;
    
    %=====================
    % Create filterbank 
    F1 = reshape(F1, Nout, nfreq);
    F2 = reshape(F2, Nout, nfreq);
    
    %=====================
    % Scalar multiply by dispersion matrix
    F1 = F1.*H;
    F2 = F2.*H;
    
    %=====================
    % If resanmpling, keep only the usable portion of the band
    if resample == 1,
        F1 = F1(floor(disc*Nout):floor((1.0-disc)*Nout));
        F2 = F2(floor(disc*Nout):floor((1.0-disc)*Nout));
    end;
    
    % Array to store de-dispersed signal
    Vc = complex(zeros(newNout, nfreq, npol, 'single'));
    
    % Inverse transform - normalise by OS factor if resample
    Vc(:,:,1) = ifft(F1, newNout);
    Vc(:,:,2) = ifft(F2, newNout);
    if resample == 1,
        Vc = Vc./os_factor;
    end;
    
    %=====================
    % Remove ends of array that are corrupted by dispersion
    z1 = Vc((1+new_n_hi):(newNout - new_n_lo), :, 1);
    z2 = Vc((1+new_n_hi):(newNout - new_n_lo), :, 2);
    
    %=====================
    % Calculate Stokes parameters: "Detection"
    sI = z1.*conj(z1) + z2.*conj(z2);
    sQ = z1.*conj(z1) - z2.*conj(z2);
    sU = 2*real(z2.*conj(z1));
    sV = 2*imag(z2.*conj(z1));
    
    % Concatenate the 4 Stokes parameters into a single matrix. This allows
    % the folding function "grpstats" to be called only once. 
    S = [sI, sQ, sU, sV];
    
    %======================
    % Fold data by calculating the running sum and running sum-of-squares
    % to each Stokes parameter assigned to a given bin
 
    % Calculate the pulsar phase for each time step. Construct a vector
    % of unique phases and compute the number of time elements assigned
    % to each phase bin. Keep a running sum of the total number of time 
    % elements in each phase bin
    
    %length(tt)
    %length(S)
    %pause
    tidx = findphase(tt, nbins, pcal);
    idx = unique(tidx);
    Vn(idx,1) = Vn(idx,1) + histc(tidx, idx);
    
    % Grouping operation to compute the sum and the sum-of-squares of
    % elements assigned to each phase bin. The functions are specified
    % as '@(x)sum(x,1)' instead of just 'sum' to enforce summation
    % over the column direction even if it only a singleton.
    [Ssum, Ssum2] = grpstats(S, tidx, {@(x)sum(x,1), @(x)sum(x.^2,1)});
    
    % Combine the sum and sum-of-squares for each Stokes parameter with 
    % the sum from previous time steps
    for s = 1:4,
        Vsum(idx,:,s) = Vsum(idx,:,s) + Ssum(:,1+(s-1)*nfreq:s*nfreq);
        Vsum2(idx,:,s) = Vsum2(idx,:,s) + Ssum2(:,1+(s-1)*nfreq:s*nfreq);
    end;
 
end;
 
fclose(fid);
 
%=================
% Package results and metadata into a data structure
 
%phase = transpose((0:(nbins-1))/(nbins-1));
phase = transpose(((1:nbins)-0.5)/nbins);
 
Vnn = repmat(Vn,1,nfreq);

Iav = Vsum(:,:,1)./Vnn;
Qav = Vsum(:,:,2)./Vnn;
Uav = Vsum(:,:,3)./Vnn;
Vav = Vsum(:,:,4)./Vnn;
 
Ierr = sqrt((Vsum2(:,:,1)./Vnn - Iav.^2)./Vnn);
Qerr = sqrt((Vsum2(:,:,2)./Vnn - Qav.^2)./Vnn);
Uerr = sqrt((Vsum2(:,:,3)./Vnn - Uav.^2)./Vnn);
Verr = sqrt((Vsum2(:,:,4)./Vnn - Vav.^2)./Vnn);
 
%Length of time covered by folded profile
Ttotal = nseries*nclip*newTout;
fprintf('Length of data analysed is %f s \n', Ttotal); 
 
dat = struct('I', Iav, 'Q', Qav, 'U', Uav, 'V', Vav, ...
             'Ierr', Ierr, 'Qerr', Qerr, 'Uerr', Uerr, 'Verr', Verr, ...
             'n', Vn, 'n_lo', new_n_lo, 'n_hi', new_n_hi, 'nclip', nclip,...
             'phase', phase, 'f', f, ...
             'nfreq', nfreq, 'nbins', nbins, 'nseries', nseries, ...
             'Tin', Tin, 'Tout', newTout, 'Nin', Nin, 'Nout', newNout, ...
             'Ttotal', Ttotal, 'pcal', pcal,...
             'df', df, 'Dconst', Dconst, 'DM', DM);
         
StokesPlots(dat);

I = mean(Iav,2);

% Plot Stokes I on a log scale
figure('name','Stokes I log scale');
set(gcf,'color','w');
plot(phase(:), 10.0*log10(I(:,1)/2.0), 'LineWidth', 1.0);
%semilogy(phase(:), Isn(:,1), 'LineWidth', 1.0);
xlabel('Phase','FontSize',12,'FontWeight','bold');
ylabel('Stokes I (dB)','FontSize',12,'FontWeight','bold');
title('Stokes I Versus Pulse Phase', 'FontSize', 14, 'FontWeight', 'bold');
axis([0 1 -Inf Inf]);
set(gca,'FontSize', 12, 'FontWeight', 'bold');


% SNR plot
Isn = zeros(nbins,1);
Ierr = mean(Ierr,2)/sqrt(nfreq);
% Calculate background level
b = phase < 0.25 | phase > 0.75; % Indices of background only phase bins
bkg = mean(I(b));
bkgerr = std(I(b))/sqrt(length(I(b)));
% Signal-to-noise
Isn(:,1) = (I - bkg)./sqrt(Ierr.^2 + bkgerr^2);
%Isn(Isn<0.0) = 0.0;
%Isn(:,1) = (I - bkg)./sqrt(Ierr.^2 + bkgerr^2 + 0.000001);

figure('name','Signal to noise');
set(gcf,'color','w');
plot(phase(:), Isn(:,1), 'LineWidth', 1.0);
%semilogy(phase(:), Isn(:,1), 'LineWidth', 1.0);
xlabel('Phase','FontSize',12,'FontWeight','bold');
ylabel('S/N','FontSize',12,'FontWeight','bold');
title('S/N Versus Pulse Phase', 'FontSize', 14, 'FontWeight', 'bold');
axis([0 1 -Inf Inf]);
%axis([phmin phmax -Inf Inf]);
set(gca,'FontSize', 12, 'FontWeight', 'bold');
 
return
end
 


function Vstream = reducebits(Vstream,Nbits)
% Mimics the effect of the data being approximated by a finite number
% of bits. This is done by extracting data that lies within an 
% optimized interval of the data (the size of the interval is given by
% off-line calculations and assumes that the maximum sigma in the data
% is 1). Extraction is achieved by scaling the data in the interval 
% to 2^32-1 integers so that data below 0 are pinned at 0 and data above
% 2^32-1 are pinned at 2^32-1. Data are converted back to floats, scaled
% and rounded to 2^Nbits-1, then scaled back to the original interval. 
 
% Optimum thresholds to minimize distortion. These are defined by separate
% simulations
bits = [2, 3, 4, 5, 6, 7, 8, 16, 32];
thres = [1.7, 2.3, 2.8, 3.2, 3.6, 4.1, 4.5, 5, 5];
a = interp1(bits,thres,Nbits,'spline');
fmin = a; % Minimum of dynamic range 
fmax = -a; % Maximum of dynamic range 
bmax = 2^Nbits-1; % Maximum integer for given number of bits
% Scale data to between 0 and 2^32-1. Values below the min are 
% pinned at min and values above the max are pinned at the max
% Convert back to float for subsequent FFT processing. 
scale32 = single(uint32((Vstream - fmin)/(fmax-fmin)*(2^32-1)));
% Scale the data to Nbits
scaleN = round(scale32/(2^32-1)*bmax);
% Rescale the data back to the original range (fmin, fmax)
Vstream = scaleN/bmax*(fmax-fmin) + fmin;
return
end
        


function [zfree, lost, SK] = RFI_remove(z0,N,M,Klo,Khi)
 
% Break the single time series z0 into M series, each with N elements
z = reshape(z0,M,N); % (M,N) matrix
 
% FFT in row direction for each of the M series. 
Z = fft(z,[],2); % (M,N) matrix
 
% Calculate instantaneous power
P = Z.*conj(Z); % (M,N) matrix
P = P./repmat(sum(P,2),1,N); %normalize each spectrum
 
% For each of N channels, calculate the sum of power over all M samples
S1 = sum(P,1); % N row vector
 
% For each of N channels, calculate the sum of squared power over all M
S2 = sum(P.^2,1); % N row vector
 
% Calculate the spectral kurtosis for each of N channels
SK = (M+1)/(M-1)*(M*S2./S1.^2 - 1); % N row vector
 
% Create logical index, with 1's identifying channels with SK's that are
% outside given limits
w = SK > Khi | SK < Klo;
 
% Record the number of channels that are outside the limit
lost = length(find(w == 1));
 
% Zero these channels across all M samples
Z(:,w) = complex(0,0); % (M,N) matrix
 
% IFFT the (M,N) matrix and reconstruct the single time series
zfree = reshape(ifft(Z,[],2),1,N*M); % N*M row vecor
 
return
end
 


function [SKlo, SKhi] = RFI_getSKlims(M,p)
% Calculates the lower (upper) limit required to make the cumulative
% probability distribution of being outside this limit be p.
% Taken from "The generalized spectral kurtosis estimator", 
% Nita & Gary (MNRAS 2010)
 
% number of power spectra that have been averaged (default to 1)
n = 1; 
 
% degree of Euler gamma function
d = 1; 
 
% Expected second moment of spectral kurtosis distribution function
m2=(2*( M^2)* (n*d) *(1 + n*d))/...
    ((-1 + M) *(6 + 5* M* (n*d) + (M^2)*( (n*d)^2)));
 
% Expected third moment of spectral kurtosis distribution function
m3=(8*( M^3)* (n*d)* (1 + n*d)* (-2 + (n*d)* (-5 + M *(4 + n*d))))/...
    (((-1 + M)^2)* (2 + M* n*d) *(3 +M* n*d)* (4 + M* n*d)* (5 + M* n*d));
 
SKlo = gammaincinv(p, 4*m2^3/m3^2, 'lower')*(m3/2/m2)+(m3-2*m2^2)/m3;
SKhi = gammaincinv(p, 4*m2^3/m3^2, 'upper')*(m3/2/m2)+(m3-2*m2^2)/m3;
 
return
end
