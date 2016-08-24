function [ pOffset,pFactor,pEdge,zEdge ] = parseGridVert( vertModel )
%PARSEGRIDVERT Retrieve vertical grid attributes
%   vertModel:      Requested model
%
%   pOffset:    Offsets for a hybrid-sigma pressure system ("Ap") [hPa]
%   pFactor:    Factors for a hybrid-sigma pressure system ("Bp") [-]
%   pEdge:      Pressure edges for a surface pressure of 1013.25 hPa [hPa]
%   zEdge:      COESA altitudes for the given pEdge [km]

% Standard GEOS-5 vertical grids
pOffsetG5 = [0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,...
             1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,...
             4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,...
             7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,...
             1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,...
             1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,...
             2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,...
             2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,...
             1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,...
             7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01,...
             4.017541e+01, 3.381001e+01, 2.836781e+01, 2.373041e+01,...
             1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,...
             9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00,...
             4.076571e+00, 3.276431e+00, 2.620211e+00, 2.084970e+00,...
             1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,...
             6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01,...
             2.113490e-01, 1.594950e-01, 1.197030e-01, 8.934502e-02,...
             6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,...
             1.000000e-02 ];
pFactorG5 = [1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,...
             9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,...
             8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,...
             7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,...
             6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,...
             4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,...
             2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,...
             6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
             0.000000e+00 ];

OKModels = {'giss40','giss23','geos5','geos5rv','era-interim','ecmwf31','ecmwf60','geos4'};
pOffset = [];
pFactor = [];
pEdge = [];
zEdge = [];
if nargin < 1
    % List available models
    fprintf('Known vertical grids:\n');
    for iModel = 1:length(OKModels)
        fprintf(' => %s\n',OKModels{iModel});
    end
    %fprintf(' => GEOS5LX (X-level interpolated GEOS-5 grid)\n');
    fprintf('Generic vertical grids:\n');
    fprintf(' => ISALEV_N_X (ISA, N levels up to X km)\n');
    fprintf(' => ISASEP_D_X (ISA, levels at separation D km up to X km)\n');
    fprintf(' => AUTO_N (automatic, N levels)\n');
    return;
elseif isempty(vertModel)
    return;
elseif length(vertModel) > 9 && strcmpi(vertModel(1:3),'ISA')
    % Fixed pressure levels
    defStr = vertModel(5:end);
    splitPt = regexp(defStr,'_');
    assert(length(splitPt) == 2,...
        'parseGridVert:badISA_1','ISA specifier cannot be parsed');
    zMax = str2double(defStr((splitPt(2)+1):end));
    switch upper(vertModel(4:6))
        case 'LEV'
            % N levels up to X km
            nLev = round(str2double(defStr((splitPt(1)+1):(splitPt(2)-1))));
        case 'SEP'
            % Separations of D km up to X km
            zDiff = str2double(defStr((splitPt(1)+1):(splitPt(2)-1)));
            nLev = round(zMax/zDiff);
            % Respecify zMax to match
            zMax = nLev * zDiff;
        otherwise
            error('parseGridVert:badISA_2','ISA specifier cannot be parsed');
    end
    zEdge = linspace(0,zMax,nLev+1)';
    pFactor = zeros(nLev+1,1);
    [~,~,pOffset] = atmoscoesa(zEdge.*1000);
    % Convert pOffset from Pa to hPa
    pOffset = pOffset.*1e-2;
    pEdge = calcPEdge(pOffset,pFactor,atmoscoesa(0)*1e-2);
    zEdge = atmospalt(pEdge.*100.0).*1.0e-3;
    return;
elseif strcmpi(vertModel,'listmodels')
    pOffset = OKModels;
    return;
elseif strncmpi(vertModel,'auto',4)
    % Auto-determine grid, if possible
    nLev = str2double(vertModel(6:end));
    % Build library of possibilities
    nModels = length(OKModels);
    gridFound = false;
    iModel = 0;
    while ~gridFound && iModel < nModels
        % Generate vertical-only grid specification
        iModel = iModel + 1;
        vertModel = OKModels{iModel};
        gSpec = genGridSpec([],vertModel);
        gridFound = gSpec.nLev == nLev;
    end
    if ~gridFound
        % Hail-mary attempt
        resMult = nLev/72;
        if abs(resMult-round(resMult)) < 1e-10
            % Enhanced GEOS-5 grid
            vertModel = ['GEOS5L',num2str(nLev)];
            [ pOffset,pFactor,pEdge,zEdge ] = parseGridVert( vertModel );
            return;
        end
    end
    assert(gridFound,'parseGridVert:noAutoMatch','Could not match a grid with %i levels to any vertical grid in library',nLev);
elseif strncmpi(vertModel,'GEOS5',5)
    % Many options
    pSfc = 1013.25;
    if strcmpi(vertModel,'GEOS5')
        % Exact match - standard vertical grid
        pOffset = pOffsetG5;
        pFactor = pFactorG5;
    elseif strcmpi(vertModel,'GEOS5RV')
        % Exact match - reduced vertical grid
        pOffset = zeros(1,48);
        pFactor = zeros(1,48);
        pOffset(1:37) = pOffsetG5(1:37);
        pFactor(1:37) = pFactorG5(1:37);
        for iLev = 37:40
            % 2 levels per RV level
            iSrc = 37 + 2*(iLev-37);
            pFactor(iLev+1) = pFactorG5(iSrc+2);
            pOffset(iLev+1) = pFactorG5(iSrc+2);
        end
        for iLev = 41:47
            % 4 levels per RV level
            iSrc = 45 + 4*(iLev-41);
            pFactor(iLev+1) = pFactorG5(iSrc+4);
            pOffset(iLev+1) = pFactorG5(iSrc+4);
        end
    elseif strncmpi(vertModel,'GEOS5L',6)
        resNum = str2double(vertModel(7:end));
        resMult = resNum/72;
        assert(abs(resMult-round(resMult))<1e-10,'parseGridVert:badResMult','Resolution does not multiply correctly');
        % Want to try and have a continuous second derivative of altitude
        % z ~ log(p)
        % Linear interpolation of layer depth (log-p)
        pFactor = ones(1,resNum+1);
        pOffset = zeros(1,resNum+1);
        
        % Logarithmic interpolation strategy (linear in altitude thickness)
        %{
        resMult = round(resMult);
        pFactor = ones(1,resNum+1);
        pOffset = zeros(1,resNum+1);
        pFactor(end) = pFactorG5(end);
        pOffset(end) = pOffsetG5(end);
        pPow = 1.0/resMult;
        hiFactor = pFactorG5(1).^pPow;
        hiOffset = pOffsetG5(1).^pPow;
        for iLev = 1:72
            loFactor = hiFactor;
            hiFactor = pFactorG5(iLev+1).^pPow;
            loOffset = hiOffset;
            hiOffset = pOffsetG5(iLev+1).^pPow;
            for subLev = 1:resMult
                outLev = (iLev-1)*resMult + subLev;
                loPow = (resMult+1-subLev);
                hiPow = subLev-1;
                newFactor = ((loFactor^loPow)*(hiFactor^hiPow));
                newOffset = ((loOffset^loPow)*(hiOffset^hiPow));
                pFactor(outLev) = newFactor;
                pOffset(outLev) = newOffset;
            end
        end
        %}
        % Linear interpolation strategy (linear in pressure)
        %{
        pFactor = interp1(0:72,pFactorG5,linspace(0,72,resNum+1));
        pOffset = interp1(0:72,pOffsetG5,linspace(0,72,resNum+1));
        %}
    else
        error('parseGridVert:badG5Model','GEOS-5 specification ''%s'' is not recognized',vertModel);
    end
    pEdge = calcPEdge(pOffset,pFactor,pSfc);
    zEdge = atmospalt(pEdge.*100.0).*1.0e-3;
    return;
elseif ~any(strcmpi(vertModel,OKModels))
    error('parseGridVert:badModel','Model ''%s'' is not recognized',...
        vertModel);
end

switch lower(vertModel)
    case {'ecmwf31'}
        % ECMWF 31-level model
        pData = [   0.000000000000e+00 0.000000000000e+00
                    2.000000000000e+01 0.000000000000e+00
                    4.000000000000e+01 0.000000000000e+00
                    6.000000000000e+01 0.000000000000e+00
                    8.000000000000e+01 0.000000000000e+00
                    9.976135361000e+01 3.910000000000e-04
                    1.182053961700e+02 2.920000000000e-03
                    1.343139392600e+02 9.194000000000e-03
                    1.473635690900e+02 2.031900000000e-02
                    1.568920745800e+02 3.697500000000e-02
                    1.626661050000e+02 5.948800000000e-02
                    1.646500573400e+02 8.789500000000e-02
                    1.629761933200e+02 1.220040000000e-01
                    1.579159860400e+02 1.614420000000e-01
                    1.498526963000e+02 2.057030000000e-01
                    1.392551785800e+02 2.541890000000e-01
                    1.266529166200e+02 3.062350000000e-01
                    1.126122887800e+02 3.611450000000e-01
                    9.771406290000e+01 4.182020000000e-01
                    8.253212096000e+01 4.766880000000e-01
                    6.761341326000e+01 5.358870000000e-01
                    5.345914240000e+01 5.950840000000e-01
                    4.050717678000e+01 6.535650000000e-01
                    2.911569385000e+01 7.105940000000e-01
                    1.954805296000e+01 7.654050000000e-01
                    1.195889791000e+01 8.171670000000e-01
                    6.381489110000e+00 8.649560000000e-01
                    2.716265450000e+00 9.077160000000e-01
                    7.206357700000e-01 9.442130000000e-01
                    0.000000000000e+00 9.729850000000e-01
                    0.000000000000e+00 9.922810000000e-01
                    0.000000000000e+00 1.000000000000e+00 ];
        pOffset = pData(end:-1:1,1);
        pFactor = pData(end:-1:1,2);
        pSfc = 1013.25;
    case {'era-interim','ecmwf60'}
        % ECMWF ERA-Interim 60-level model
        pData = [   0.000000000000e+00 1.000000000000e+00
                    0.000000000000e+00 9.976300000000e-01
                    7.367740000000e+00 9.940190000000e-01
                    6.588920000000e+01 9.882700000000e-01
                    2.103940000000e+02 9.796630000000e-01
                    4.673330000000e+02 9.676450000000e-01
                    8.553620000000e+02 9.518220000000e-01
                    1.385910000000e+03 9.319400000000e-01
                    2.063780000000e+03 9.078840000000e-01
                    2.887700000000e+03 8.796570000000e-01
                    3.850910000000e+03 8.473750000000e-01
                    4.941780000000e+03 8.112530000000e-01
                    6.144320000000e+03 7.715970000000e-01
                    7.438800000000e+03 7.287860000000e-01
                    8.802360000000e+03 6.832690000000e-01
                    1.020950000000e+04 6.355470000000e-01
                    1.163280000000e+04 5.861680000000e-01
                    1.304320000000e+04 5.357100000000e-01
                    1.441110000000e+04 4.847720000000e-01
                    1.570640000000e+04 4.339630000000e-01
                    1.689950000000e+04 3.838920000000e-01
                    1.796140000000e+04 3.351550000000e-01
                    1.886480000000e+04 2.883230000000e-01
                    1.958430000000e+04 2.439330000000e-01
                    2.009740000000e+04 2.024760000000e-01
                    2.038450000000e+04 1.643840000000e-01
                    2.042990000000e+04 1.300230000000e-01
                    2.022220000000e+04 9.967470000000e-02
                    1.975510000000e+04 7.353380000000e-02
                    1.902770000000e+04 5.169040000000e-02
                    1.804520000000e+04 3.412120000000e-02
                    1.681950000000e+04 2.067790000000e-02
                    1.537980000000e+04 1.114290000000e-02
                    1.377530000000e+04 5.081120000000e-03
                    1.207740000000e+04 1.815160000000e-03
                    1.037612000000e+04 4.613950000000e-04
                    8.765050000000e+03 7.582350000000e-05
                    7.306630000000e+03 0.000000000000e+00
                    6.018020000000e+03 0.000000000000e+00
                    4.906710000000e+03 0.000000000000e+00
                    3.960290000000e+03 0.000000000000e+00
                    3.196420000000e+03 0.000000000000e+00
                    2.579890000000e+03 0.000000000000e+00
                    2.082270000000e+03 0.000000000000e+00
                    1.680640000000e+03 0.000000000000e+00
                    1.356470000000e+03 0.000000000000e+00
                    1.094830000000e+03 0.000000000000e+00
                    8.836600000000e+02 0.000000000000e+00
                    7.132180000000e+02 0.000000000000e+00
                    5.756510000000e+02 0.000000000000e+00
                    4.646180000000e+02 0.000000000000e+00
                    3.739720000000e+02 0.000000000000e+00
                    2.984960000000e+02 0.000000000000e+00
                    2.347790000000e+02 0.000000000000e+00
                    1.805840000000e+02 0.000000000000e+00
                    1.344830000000e+02 0.000000000000e+00
                    9.563700000000e+01 0.000000000000e+00
                    6.364780000000e+01 0.000000000000e+00
                    3.842530000000e+01 0.000000000000e+00
                    2.000000000000e+01 0.000000000000e+00
                    0.000000000000e+00 0.000000000000e+00];
        pOffset = pData(:,1).*0.01;
        pFactor = pData(:,2);
        pSfc = 1013.25;
    case {'giss40'}
        % Define SIGMA edges from GISS 40L model for GEOS-Chem compatability
        nLev = 40;
        sigEdge = [ 1.00000000000e0, 0.97601918465e0, 0.94964028777e0,...
                    0.91966426859e0, 0.88729016787e0, 0.85131894484e0,...
                    0.80935251799e0, 0.76139088729e0, 0.70743405276e0,...
                    0.64988009592e0, 0.58992805755e0, 0.52877697842e0,...
                    0.46642685851e0, 0.40647482014e0, 0.34892086331e0,...
                    0.29496402878e0, 0.24460431655e0, 0.19904076739e0,...
                    0.15827338129e0, 0.12110311751e0, 0.08752997602e0,...
                    0.05635491607e0, 0.02757793765e0, 0.00000000000e0,...
                   -0.02637889688e0,-0.05035971223e0,-0.07194244604e0,...
                   -0.09232613909e0,-0.11151079137e0,-0.12829736211e0,...
                   -0.14268585132e0,-0.15587529976e0,-0.16786570743e0,...
                   -0.17311750600e0,-0.17606714628e0,-0.17772182254e0,...
                   -0.17865707434e0,-0.17918225420e0,-0.17947721823e0,...
                   -0.17964268585e0,-0.17973621103e0                   ];
        
        lTransition = 24;
        pTransition = 150.0;
        pSfcFix = 984;
        pOffset = zeros(nLev,1);
        pFactor = zeros(nLev,1);
        for iLev = 1:(nLev+1)
            if iLev >= lTransition
                pOffset(iLev) = sigEdge(iLev)*(pSfcFix-pTransition) + pTransition;
                pFactor(iLev) = 0.0;
            else
                pOffset(iLev) = pTransition*(1.0-sigEdge(iLev));
                pFactor(iLev) = sigEdge(iLev);
            end
        end
        
        pSfc = 1013.25;
    case {'giss23'}
        % Define SIGMA edges from GISS 23L model for GEOS-Chem compatability
        nLev = 23;
        sigEdge =   [1.0e0,           0.9712230e0,     0.9340528e0,...
                     0.8800959e0,     0.8021583e0,     0.6714628e0,...
                     0.5035971403e0,  0.3297362030e0,  0.1966426820e0,...
                     0.1139088720e0,  0.0503597111e0,  0.0000000000e0,...
                    -0.0395683460e0, -0.0764988065e0, -0.1124700233e0,...
                    -0.1419664323e0, -0.1585131884e0, -0.1678657085e0,...
                    -0.1743045598e0, -0.1781055182e0, -0.1793033630e0,...
                    -0.1796822548e0, -0.1798187047e0, -0.1798536479e0 ];
        
        lTransition = 12;
        pTransition = 150.0;
        pSfcFix = 984;
        pOffset = zeros(nLev,1);
        pFactor = zeros(nLev,1);
        for iLev = 1:(nLev+1)
            if iLev >= lTransition
                pOffset(iLev) = sigEdge(iLev)*(pSfcFix-pTransition) + pTransition;
                pFactor(iLev) = 0.0;
            else
                pOffset(iLev) = pTransition*(1.0-sigEdge(iLev));
                pFactor(iLev) = sigEdge(iLev);
            end
        end
        
        pSfc = 1013.25;
    case {'geos5rv'}
        pSfc = 1013.25;
        pOffset = [  0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,...
                     1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,...
                     4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,...
                     7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,...
                     1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,...
                     1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,...
                     2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,...
                     2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,...
                     1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,...
                     7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01,...
                     1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00,...
                     6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02 ];
        pFactor = [  1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,...
                     9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,...
                     8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,...
                     7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,...
                     6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,...
                     4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,...
                     2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,...
                     6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,...
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,...
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00 ];
    case {'geos5'}
        pSfc = 1013.25;
        pOffset = pOffsetG5;
        pFactor = pFactorG5;
    case {'geos4'}
        pOffset = [  0.000000e+0;  0.000000e+0;...
                    12.704939e+0; 35.465965e+0;...
                    66.098427e+0;101.671654e+0;...
                   138.744400e+0;173.403183e+0;...
                   198.737839e+0;215.417526e+0;...
                   223.884689e+0;224.362869e+0;...
                   216.864929e+0;201.192093e+0;...
                   176.929993e+0;150.393005e+0;...
                   127.837006e+0;108.663429e+0;...
                    92.365662e+0; 78.512299e+0;...
                    66.603378e+0; 56.387939e+0;...
                    47.643932e+0; 40.175419e+0;...
                    33.809956e+0; 28.367815e+0;...
                    23.730362e+0; 19.791553e+0;...
                    16.457071e+0; 13.643393e+0;...
                    11.276889e+0;  9.292943e+0;...
                     7.619839e+0;  6.216800e+0;...
                     5.046805e+0;  4.076567e+0;...
                     3.276433e+0;  2.620212e+0;...
                     2.084972e+0;  1.650792e+0;...
                     1.300508e+0;  1.019442e+0;...
                     0.795134e+0;  0.616779e+0;...
                     0.475806e+0;  0.365041e+0;...
                     0.278526e+0;  0.211349e+0;...
                     0.159495e+0;  0.119703e+0;...
                     0.089345e+0;  0.066000e+0;...
                     0.047585e+0;  0.032700e+0;...
                     0.020000e+0;  0.010000e+0];
        pFactor  = [ 1.000000e+0; 0.985110e+0; ...
            0.943290e+0; 0.867830e+0; ...
            0.764920e+0; 0.642710e+0; ...
            0.510460e+0; 0.378440e+0; ...
            0.270330e+0; 0.183300e+0; ...
            0.115030e+0; 0.063720e+0; ...
            0.028010e+0; 0.006960e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0; ...
            0.000000e+0; 0.000000e+0 ];
        pSfc = 1013.25;

    otherwise
        % No vertical grid association
        pEdge = [];
        pFactor = [];
        zEdge = [];
        pOffset = [];
end

if ~isempty(pFactor)
    pEdge = calcPEdge(pOffset,pFactor,pSfc);
    zEdge = atmospalt(pEdge.*100.0).*1.0e-3;
end

end

