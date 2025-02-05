%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AE302 Lab 1 - Wyatt Welch
%
% NOTE: For groups 1 and 2, there are seperate and 
% nearly identical documents. This document has the 
% resulting data sets from those two imported for 
% simplicity and readability.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all 

% Raw Data
barIn = [29.2 29.2 29.3 29.4 29.4 29.3 29.3 29.8 29.9 29.1 29.4 29.2 29.4 29.9 29.1 29.9 29.4 29.7 29.6 29.2];
barMm = [760.5 760.6 760.7 760.3 761.1 760.4 760.3 760.6 760.0 760.5 761.2 760.0 760.4 760.1 760.3 760.1 760.6 759.2 758.0 759.8];
temp = [73.0 73.0 74.0 74.0 74.0 73.0 73.0 73.0 73.0 73.0 75.0 74.0 73.0 73.0 73.0 73.0 73.0 74.0 74.0 73.0];
tempV = temp';
rawtable3 = [barIn;barMm;temp]';

WbarIn = 29.4;
WbarMm = 760.6;
Wtemp = 73;

% Constants
gm = 9.8066; %m/s^2
gf = 32.174; %ft/s^2
R = 8.314; %J/molK
M = .028964; %kg/mol


t = 2.093; %Using table

% Latitude and Temperature Corrections
tAdjust = [.117 .117 .120 .120 .120 .118 .118 .120 .120 .121 .123 .119 .118 .121 .117 .121 .118 .121 .122 .117];
tCor = temp - tAdjust;
tCorK = (tCor - 32) .* (5/9) + 273.15;

WtCor = 72.882;
WtCorC = 22.1722;
WtCorK = 295.862;

pAdjust = [.032 .032 .032 .032 .032 .032 .032 .033 .033 .032 .032 .032 .032 .033 .032 .033 .032 .033 .033 .032];
pCor = barIn - pAdjust;
pCorPa = pCor .* 3386.39;
WpCor = 29.3682;
WpCorAtm = .98152;
WpCorPa = 99452.5;
combtp3 = [pCorPa;tCorK]';

%%%%%%%% Calculations %%%%%%%%
N = length(temp);

tSampleMean = sum(tCorK)./N;
pSampleMean = sum(pCorPa)./N;

tSampleSD = sqrt(sum((tCorK-tSampleMean).^2) ./ (N-1));
tSampleSDpm = tSampleSD * t;
pSampleSD = sqrt(sum((pCorPa-pSampleMean).^2) ./ (N-1));
pSampleSDpm = pSampleSD * t;

tSampleSDMean = tSampleSD ./ sqrt(N);
pSampleSDMean = pSampleSD ./ sqrt(N);

% Calcuate Density
Wdensity = (WpCorPa * M) / (R * WtCorK);
density = (pCorPa .* M) ./ (R .* tCorK);

% Calculate Confidence Intervals
tCI = [tSampleMean + (t*tSampleSD)/sqrt(N), tSampleMean - (t*tSampleSD)/sqrt(N)];
pCI = [pSampleMean + (t*pSampleSD)/sqrt(N), pSampleMean - (t*pSampleSD)/sqrt(N)];

wtCI = [WtCorK + (t*tSampleSD)/sqrt(N), WtCorK - (t*tSampleSD)/sqrt(N)];
wpCI = [WpCorPa + (t*pSampleSD)/sqrt(N), WpCorPa - (t*pSampleSD)/sqrt(N)];

dMean = sum(density) ./ N;
dSD = sqrt(sum((density-dMean).^2) ./ (N-1));
dSDMean = dSD ./ sqrt(N);
dCI = [dMean + (t*dSD)/sqrt(N), dMean - (t*dSD)/sqrt(N)];
wdCI = [Wdensity + (t*dSD)/sqrt(N), Wdensity - (t*dSD)/sqrt(N)];

% Running Averages
running_avg = zeros(1, N);

cumulative_sum = 0; 
for i = 1:N
    cumulative_sum = cumulative_sum + pCorPa(i);
    running_avg(i) = cumulative_sum / i;
end

std_devs = zeros(1, N); % Preallocate array for standard deviations

for i = 1:N
    std_devs(i) = sqrt(sum((pCorPa(1:i)-pSampleMean).^2) ./ (N-1)); % Sample standard deviation (normalizing by N-1)
end


sd_means = zeros(1, N); % Preallocate array for SD/sqrt(N)

for i = 1:N
    sd_means(i) = sqrt(sum((pCorPa(1:i)-pSampleMean).^2) ./ (N-1)) / sqrt(i); % Standard deviation of the means
end

%%%%%%%% Plotting %%%%%%%%
running_avg1 = [101212.424320000	101263.220170000	101280.152120000	101305.550045000	101314.016020000	101314.016020000	101299.502920000	101288.194796250	101294.450211111	101293.020402000	101298.007630909	101299.341663333	101300.470460000	101306.033815000	101308.823555333	101308.936435000	101311.227228235	101316.838011667	101320.075875790	101311.306908000]./1000;
running_avg2 = [101280.152120000	101263.220170000	101246.288220000	101246.288220000	101212.424320000	101229.356270000	101231.775120000	101212.001021250	101211.671788889	101201.587872000	101211.808612727	101211.859921667	101217.113167692	101216.536365000	101211.747042000	101215.810710000	101223.579487059	101230.296933889	101232.742660000	101228.340353000]./1000;
running_avg3 = running_avg./1000;

figure(1)
plot(1:N, running_avg1,'-o')
grid on, hold on
plot(1:N, running_avg2,'-square')
plot(1:N, running_avg3,'-^')
xlabel("Number of Samples")
ylabel("Sample Mean (kPa)")
title("Number of Samples vs. Sample mean")
legend("Group 1",'Group 2','Group 3')


std1 = [22.6852214930422	22.6937337715754	22.7022428584028	27.8660360685761	29.1018072906668	29.1084432062223	36.9042103299644	43.7308983398862	44.3887075507091	44.9604566713035	45.7366580053691	45.7408806627818	45.7451029304079	48.2621388606169	48.9860517781275	48.9862981986143	49.6996643108532	54.8273829910075	56.9443722498587	68.5832534711971];
std2 = [11.8864345494337	12.5794007252270	13.0986243471517	13.7305482809560	37.3429359506760	42.1998518074240	42.4002534058043	55.3019910041914	55.4790038545627	61.6819280963565	64.7378782540512	64.8407704921704	65.9212624599702	66.0698305127477	68.8000605777787	69.6912497385152	74.8929915879104	79.4924041779894	80.2649629452783	82.5269053589780];
std3 = std_devs;
figure(2)
plot(1:N, std1,'-o')
grid on, hold on
plot(1:N, std2,'-square')
plot(1:N, std3,'-^')
xlabel("Number of Samples")
ylabel("Sample Standard Deviation")
title("Number of Samples vs. Sample Standard Deviation")
legend("Group 1",'Group 2','Group 3')


stdm1 = [22.6852214930422	16.0468930403231	13.1071460255071	13.9330180342880	13.0147238740060	11.8834721770047	13.9484804091867	15.4612073817565	14.7962358502364	14.2177447722632	13.7901212516922	13.2042548818138	12.6874087859203	12.8985991692383	12.6481441822095	12.2465745496536	12.0539391477287	12.9229381025511	13.0639349494998	15.3356816879695];
stdm2 = [11.8864345494337	8.89497955607097	7.56249429284182	6.86527414047802	16.7002686530264	17.2280173582092	16.0257894339825	19.5522064260906	18.4930012848542	19.5055383255221	19.5192047174905	18.7179181490586	18.2832686109802	17.6579049557803	17.7640992557199	17.4228124346288	18.1642185256147	18.7365393490260	18.4140453781845	18.4535770355366];
stdm3 = sd_means;
figure(3)
plot(1:N, stdm1,'-o')
grid on, hold on
plot(1:N, stdm2,'-square')
plot(1:N, stdm3,'-^')
xlabel("Number of Samples")
ylabel("Sample Standard Deviation of the Means")
title("Number of Samples vs. Sample Standard Deviation of the Means")
legend("Group 1",'Group 2','Group 3')