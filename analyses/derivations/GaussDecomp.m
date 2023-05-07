function [amps, compspec, sumspec] = GaussDecomp(wl,ap,wlunc,apunc)

% AC-S "un-smoothing" and spectral decomposition method
%
% Ron Zaneveld, WET Labs, Inc., 2005
% Ali Chase, University of Maine, 2014
% 
% Vectorized: Guillaume Bourdin, University of Maine, 2020
%
% See the following publication for details on the method:
% Chase, A., et al., Decomposition of in situ particulate absorption
% spectra. Methods in Oceanography (2014), http://dx.doi.org/10.1016/j.mio.2014.02.022
%
% INPUTS:
% wl     -  wavelengths associated with measured particulate absorption
%           spectra
% ap     -  measured particulate absorption spectra
% wlunc  -  wavelengths associated with uncertainty of measured
%           particulate absorption spectra
% apunc  -  uncertainty of measured particulate absorption spectra
% acs    -  1 (for AC-S data) or 0 (for a different instrument). This code
%           was designed for use with particulate absorption spectra measured
%           using a WETLabs AC-S instrument, and thus it applies a spectral
%           un-smoothing that is specific to that instrument. To use with
%           other types of absorption data, input "0" and the correction
%           applied for AC-S will be bypassed.
%
% The uncertainty values we use are the standard deviation of one-minute
% binned particulate absorption spectra. If the uncertainty is unknown,
% then all wavelenghts will be treated equally, i.e. no spectral weighting
% will be applied. In this case populate apunc with -999 values.
%
% OUTPUTS:
% amps      -  the amplitudes of 12 Gaussian function peaks and one non-algal 
%              particle (NAP) function
% compspec  -  the component spectra after spectral decomposition; the
%              Gaussian and NAP functions multiplied by the 'amps'
%              values
% sumspec   -  the sum of the component spectra. will be simliar to the
%              original measured spectrum, but the fit may not be as good in
%              parts of the spectrum with higher uncertainty.
%
% The amplitudes of the Gaussian functions represent the relative amounts
% of different phytoplankton pigments and can be compared to HPLC pigment
% concentrations to evaluate the method. See Chase et al. (2014) for
% detail on such evaluation (reference above).


% Peak center values ("peak_loc") determined using a interative
% approach that allows the location to vary (uses the matlab
% function LSQNONLIN), and are rounded to nearest integer.
% Sigma values ("lsqsig") are determined similarly. FWHM = sigma*2.355
abscorr = ap;

peak_loc=[406,434,453,470,492,523,550,584,617,638,660,675];

lsqsig=[16,12,12,13,16,14,14,16,13,11,11,10];

onenm = 400:1:720;

fprintf('Gaussian decomposition ... ')

%interpolate the un-smoothed ap spectra to one nm resolution
acorr2onenm = interp1(wl, abscorr', onenm, 'spline');

%define the matrix of component Gaussian functions using the peaks
%and widths (sigma) above
coef2 = exp(-0.5 .* (((onenm .* ones(size(peak_loc,2),1))' - peak_loc .* ...
    ones(size(onenm,2),1)) ./ lsqsig) .^ 2);

%define a function for non-algal particles and concatenate this to the Gaussian matrix
coef2nap = exp(-0.01 * (onenm - 400));
coef2 = [coef2nap', coef2];

%normalize both the component functions and the measured ap
%spectrum by the uncertainty (std dev) in the ap spectra
apunc_int = interp1(wlunc, apunc', onenm, 'linear', 'extrap');
acorr2onenm_new = (acorr2onenm ./ apunc_int);

amps = NaN(size(ap,1), size(coef2,2));
sumspec_temp = NaN(size(coef2,1), size(ap,1));
compspec_temp = NaN(size(coef2,1), size(coef2,2), size(ap,1));

for i = 1:size(acorr2onenm_new,2)
    coef2_new = coef2 ./ apunc_int(:,i);
    %Inversion analysis
    amps(i, :) = lsqnonneg(coef2_new, acorr2onenm_new(:, i));
    
	% Using the inverted amplitudes, build a matrix of new component
    % spectra (compspec) and the sum of the Gaussian and nap functions (sumspec):
    sumspec_temp(:, i) = sum(amps(i, :) .* coef2, 2);
    compspec_temp(:, :, i) = amps(i, :) .* coef2;
end

%interpolate back to the original resolution
compspec = interp1(onenm', compspec_temp, wl, 'spline');
sumspec = interp1(onenm, sumspec_temp, wl, 'spline')';
fprintf('Done\n')
end


