function [ap_corr, cp_corr] = ResidualTemperatureAndScatteringCorrection(ap, cp, wl)
% Function from Emmanuel Boss improved by Nils Haëntjens

% Sullivan et al. 2006 values
psi_wl = [400;402;404;406;408;410;412;414;416;418;420;422;424;426;428;430;432;434;436;438;440;442;444;446;448;450;452;454;456;458;460;462;464;466;468;470;472;474;476;478;480;482;484;486;488;490;492;494;496;498;500;502;504;506;508;510;512;514;516;518;520;522;524;526;528;530;532;534;536;538;540;542;544;546;548;550;552;554;556;558;560;562;564;566;568;570;572;574;576;578;580;582;584;586;588;590;592;594;596;598;600;602;604;606;608;610;612;614;616;618;620;622;624;626;628;630;632;634;636;638;640;642;644;646;648;650;652;654;656;658;660;662;664;666;668;670;672;674;676;678;680;682;684;686;688;690;692;694;696;698;700;702;704;706;708;710;712;714;716;718;720;722;724;726;728;730;732;734;736;738;740;742;744;746;748;750]';
psiT = [0.0001;0.0001;0.0001;0.0001;0;0;0;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0003;0.0003;0.0004;0.0005;0.0006;0.0006;0.0007;0.0008;0.0009;0.001;0.001;0.001;0.001;0.001;0.0009;0.0009;0.0008;0.0007;0.0006;0.0006;0.0005;0.0004;0.0003;0.0003;0.0002;0.0001;0.0001;0;0;0;0;0;0;0;0;0;0.0001;0.0001;0.0001;0.0002;0.0002;0.0002;0.0001;0.0001;0.0001;0;0;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;-0.0001;0;0;0.0001;0.0002;0.0003;0.0005;0.0007;0.0009;0.0013;0.0017;0.0021;0.0026;0.0032;0.0038;0.0045;0.0054;0.0063;0.0073;0.0083;0.0094;0.0104;0.0113;0.0121;0.0128;0.0133;0.0136;0.0136;0.0133;0.0129;0.0124;0.0116;0.0107]';
psiT = interp1(psi_wl, psiT, wl);

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
% Find nearest wavelength to greater than 715 nm to use as reference for correction
iref = find(715 <= wl, 1,'first'); % 730
% If ACS spectrum does not go up to 730 nm take the closest wavelength to 730 nm
if isempty(iref); [~, iref] = max(wl); end % works as there is data in iNIR so lowest wavelength is 710


% Initialize output arrays
deltaT = NaN(size(ap,1),1);

% Init routine parameters
bp = cp - ap;

% Run minimization routine on good spectrum only
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

