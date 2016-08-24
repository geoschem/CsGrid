function [T_K, a_ms, P_Pa, rho_kgm3] = atmosisa(z_m)
%ATMOSISA Lightweight version of atmosisa from the Aerospace Toolbox

% COESA data
HVec = [-0.1,0,11,20,32,47,51,71,84.8520,1e6]; % km
LVec = [0.0,-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0,0.0]; % K/km
TRef = 288.15;
HDiff = diff(HVec);

% Change in temperature between each level
TDiff = HDiff.*LVec;
TVec = cumsum([TRef,TDiff]);

% Gas composition
gasMW = [	28.0134,...
			31.9988,...
			39.948,...
			44.00995,...
			20.183,...
			4.0026,...
			83.80,...
			131.30,...
			16.04303,...
			2.01594];

gasFrac = [	0.78084,...
			0.209476,...
			0.00934,...
			0.000314,...
			0.00001818,...
			0.00000524,...
			0.00000114,...
			0.000000087,...
			0.000002,...
			0.0000005];

% Normalize, to be 100% safe
gasFrac = gasFrac./sum(gasFrac);

T_K = zeros(size(z_m));
a_ms = zeros(size(z_m));
P_Pa = zeros(size(z_m));
rho_kgm3 = zeros(size(z_m));

T_K = interp1(HVec,TVec,z_m./1000.0);

RStar = 8.31432;
NAvog = 6.022169e23; % molec/mol

%dynVisc = (T_K.^1.5 .* 1.458e-6)./(T_K + 110.4);
%kinVisc = dynVisc./rho_kgm3;

g0 = 9.80665; %m/s2

iLo = 1;
iHi = 2;
zLo = HVec(1) .* 1000;
zHi = HVec(2) .* 1000;
MgR = sum(gasFrac.*gasMW).*1e-3.*g0./RStar;
TLo = TVec(iLo);
alphaTemp = 0;
% Exponential offset
PBase = 101325 .*exp(-MgR.*zLo./TLo);
for iPoint = 1:length(T_K)
	zCurr = z_m(iPoint);
	while zCurr > zHi
		%disp([TVec(iHi) TVec(iLo) MgR alphaTemp zLo zHi PBase])
		if abs(alphaTemp) > 0
			PNew = PBase .* (TVec(iHi)./TVec(iLo)).^(MgR./-alphaTemp);
		else
			PNew = PBase .* exp(MgR.*(zLo-zHi)./TLo);
		end
		%fprintf('%5.2f km, %5.2f K, %9.2f -> %9.2f hPa, %8.5f K/m,\n',HVec(iLo),TVec(iLo),PBase./100,PNew./100,alphaTemp);
		PBase = PNew;
		iLo = iHi;
		iHi = iHi + 1;
		zLo = zHi;
		zHi = HVec(iHi) .* 1000;
		TLo = TVec(iLo);
		alphaTemp = LVec(iLo)./1000;
	end
	if abs(alphaTemp) > 0
		P_Pa(iPoint) = PBase.*((T_K(iPoint)./TLo)).^(MgR./-alphaTemp);
	else
		P_Pa(iPoint) = PBase.*exp(MgR.*(zLo-zCurr)./TLo);
	end
end

end
