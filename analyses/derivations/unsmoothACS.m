function acs_unsmoothed = unsmoothACS(spec, lambda)
% AC-S "un-smoothing" and spectral decomposition method
% Ron Zaneveld, WET Labs, Inc., 2005
% Ali Chase, University of Maine, 2014
% Modified by Guillaume Bourdin, University of Maine, 2020
%
% Unsmoothing method from spectral decomposition in:
% Chase, A., et al., Decomposition of in situ particulate absorption
% spectra. Methods in Oceanography (2014), http://dx.doi.org/10.1016/j.mio.2014.02.022
%
% Input:
%    - spec: <NxM double> ACs spectra, either absorption or attenuation
%    - lambda: <1xM double> wavelengths
%
% Output:
%    - acs_unsmoothed: <NxM double> ACs spectra unsmoothed
%
% WARNING: avoid NaNs in spectra to prevent unrealistic interpolation and extrapolation
%%
fprintf('Unsmoothing ... ')

% force lambda in columns
lambda = lambda(:)';

% Set up filter factors at every 0.1 nm from 1 to 799 nm, with center
% wavelength at centwavel (i.e. at the data wavelengths)
wavelength = .1:.1:799; % Thus index of 1 nm = 10; 356 nm= 3560;
SIG1 = (-9.845*10^-8 .* lambda.^3 + 1.639*10^-4 * lambda.^2 - 7.849*10^-2 * lambda + 25.24) / 2.3547 ;
for j = 1:size(lambda,2)
  for jkl = 1:size(wavelength,2)
    filtfunc(jkl,j) = (1/(sqrt(2*pi)*SIG1(j)))*exp(-0.5*((wavelength(jkl)-lambda(j))/SIG1(j)).^2); % First term normalizes area under the curve to 1.
  end
end

% Convolve the measurement with the fiter factors add the difference to
% the measured spectrum to get the first corrected spectrum.
% This is the corrected absorption spectrum "ap".
minwavel = min(lambda);
maxwavel = max(lambda);

centwavel = minwavel:.1:maxwavel;% The range of centwavel is 0.1 nm.
splinap = spline(lambda, spec, centwavel); % Spline the measured data to every 0.1 nm.
% We need data from 0 to 799 nm to multiply by filtfac.
absspec = zeros(size(spec,1), size(wavelength,2));
absspec(:, minwavel*10:maxwavel*10) = splinap;
absspec(:, 1:minwavel*10-1) = ones(1, size(1:minwavel*10-1,2)) .* absspec(:, minwavel*10);
aspecprime = absspec';

meassignal6 = NaN(size(aspecprime, 2), size(lambda, 2));
for j = 1:size(aspecprime, 2)        
  measur2 = aspecprime(:,j) .* filtfunc; % the measured signal for every filter factor.
  meassignal6(j,:) = 0.1 * sum(measur2); % The measured spectrum at a wavelength i is the sum of what a filter measured at
end
acs_unsmoothed = spec - meassignal6 + spec;
fprintf('Done\n')
end