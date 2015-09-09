function phase = findphase(t, nbins, pcal)
% Determines the phase from the absolute time, 't', the number of bins
% per pulse period 'n', and the pulsar calibration structure 'pcal'.
% 
% Inputs:
% --------
%   t       - vector of times (in seconds)
%   nbins   - number of bins per pulsar period
%   pcal    - structure containing pulsar details
% 
% Outputs:
% --------
%
%    dat -  pulsar "phase", given as integer between 1 and nbins
%
% Description:
% ------------
% Calculates the pulsar phase corresponding to the requested times 
% by using a known pulsar calibration
% 
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% D. Hicks         21-Apr-2014  Original version
% ----------------------------------------------------------------------
 
%phase = 1 + transpose(round(mod(t, pcal.a)/pcal.a*(nbins-1)));
 
% fractional period
f = mod(t, pcal.a)/pcal.a; 
 
% bin assignments are the center of each bin (Matlab indexing).
phase = transpose(round(f*nbins+0.5));
 
 
return
end
