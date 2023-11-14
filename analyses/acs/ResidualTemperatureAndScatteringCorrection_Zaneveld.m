%% Residual Temperature And Scattering Correction (Zaneveld)
function [ap_corr, cp_corr] = ResidualTemperatureAndScatteringCorrection_Zaneveld(ap, cp, wl, psi)
  % Function from Emmanuel Boss improved by Nils Haëntjens
  % Assumes negligible ap in NIR after Zaneveld 1994 method 3
  psiT = interp1(psi.wl, psi.psiT, wl);
  
  % Parameters of minization routine
  opts = optimset('fminsearch');
  % opts = optimset(opts,'NonlEqnAlgorithm', 'gn'); % Does not work on R2017a
  opts = optimset(opts,'MaxIter',20000000);
  opts = optimset(opts,'MaxFunEvals',20000);
  opts = optimset(opts,'TolX',1e-8);
  opts = optimset(opts,'TolFun',1e-8);
  
  % Find Near Infrared & references
  iNIR = 710 <= wl &  wl <= 750;  % spectral srange for optimization (710 to 750nm)
  if isempty(iNIR); error('Unable to perform correction as no wavelength available in NIR.'); end
  % Find nearest wavelength to greater than 730 nm to use as reference for correction
  iref = find(730 <= wl, 1,'first'); % 715 730
  % If ACS spectra do not go up to 730 nm take the closest wavelength to 730 nm
  if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710
  
  % Initialize output arrays
  deltaT = NaN(size(ap,1),1);
  
  % Init routine parameters
  bp = cp - ap;
  
  % Run minimization routine on good spectra only
  sel = find(all(isfinite(ap),2));
  for k = sel'
    deltaT(k) = fminsearch(@costFun_RTSC, 0, opts, ap(k,:), bp(k,:), psiT, iNIR, iref);         
  end
  ap_corr = ap - psiT.*deltaT - ((ap(:,iref) - psiT(iref).*deltaT) ./ bp(:,iref)) .* bp; 
  cp_corr = cp - psiT.*deltaT;
end

function cost = costFun_RTSC(deltaT, ap, bp, psiT, iNIR, iref)
    cost = sum(abs(ap(iNIR) - psiT(iNIR).*deltaT - ((ap(iref)-psiT(iref).*deltaT)./bp(iref)).*bp(iNIR)));
end

