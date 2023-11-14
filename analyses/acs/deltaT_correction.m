function [a_Tcorr, c_Tcorr] = deltaT_correction(a, c, wl, psi)

  opts = optimset('fminsearch');      
  opts = optimset(opts,'MaxIter',20000000); 
  opts = optimset(opts,'MaxFunEvals',20000);
  opts = optimset(opts,'TolX',1e-8);
  opts = optimset(opts,'TolFun',1e-8);
  
  phi_T=interp1(psi.wl, psi.psiT, wl, 'linear');
  a_Tcorr = NaN(size(a));
  c_Tcorr = NaN(size(c));
  for i = 1:size(a, 1)
    % correct a
    spectra = a(i,:);
    I = find(wl>=710 & wl<=750);  % spectral range for optimization (710 to 750nm)   
    delT = 0;
    offset = 0;
    [x1] = fminsearch(@f_T, [delT, offset], opts, spectra(I), phi_T(I)');
    a_Tcorr(i,:)=spectra-x1(1).*phi_T';
    % correct c
    spectra = c(i,:);
    delT = 0;
    offset = 0;
    [x1] = fminsearch(@f_T, [delT, offset], opts, spectra(I), phi_T(I)');
    c_Tcorr(i,:) = spectra-x1(1).*phi_T';
  end
end

function costf = f_T(x0, spectra, psiT)
  costf = sum((spectra - psiT.*x0(1) -x0(2)).^2);
end