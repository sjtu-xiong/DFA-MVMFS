function [leaders,scales,ncount] = dwtleaders(x,Lo,Hi,Ncfs)
% This function is for internal use only and may change in a future
% release.
%
% The wavelet leader algorithm used here is due to Dr. Stephane Roux
% and colleagues:
% http://www.ens-lyon.fr/PHYSIQUE/teams/signaux-systemes-physique
% http://perso.ens-lyon.fr/stephane.roux/
%
% and is described in these references:
%
% Herwig Wendt, Patrice Abry, 
% "Multifractality Tests using Bootstrapped Wavelet Leaders", 
% IEEE Trans. Signal Processing, vol. 55, no. 10, pp. 4811-4820, 2007.
%
% Wendt,H., Abry, P. & Jaffard, S. (2007)
% "Bootstrap for Empirical Multifractal Analysis", IEEE Signal Processing
% Magazine, 24, 4, 38-48. 
%
% Wendt,H (2008) "Contributions of Wavelet Leaders and 
% Bootstrap to Multifractal Analysis: Images, Estimation Performance,
% Dependence Structure and Vanishing Moments. Confidence Intervals 
% and Hypothesis Tests."

%   Copyright 2016-2022 The MathWorks, Inc.
%#codegen

% Obtain the length of the scaling filter
nwav = coder.const(numel(Lo));
nvalid = coder.const(nwav-1);
% Filter should not be longer than 300 taps for codegen
if ~isempty(coder.target) && ~strcmp(eml_option('UseMalloc'),'VariableSizeArrays')
    coder.internal.assert(nwav <= 300,'Wavelet:FunctionInput:FilterLenDMAOff');
end
% Obtain the length of the signal and determine the level of the DWT
n = Ncfs-nvalid;
Nlevels = min(fix(log2(n/(nwav+1))),fix(log2(n)));
% This makes sure we have at least three levels but due to boundary effects
% we still may not have enough leaders at the end, so we will check again
% before returning
coder.internal.errorIf(Nlevels < 3, 'Wavelet:mfa:InsufficientLeaders');
% Shifts for wavelet and scaling coefficients
x0=2;
if ~isempty(coder.target) && ...
        (coder.internal.isConst(Ncfs) || ...
        ~coder.internal.eml_option_eq('UseMalloc','VariableSizeArrays'))
    UB = coder.const(Ncfs);
else
    UB = coder.const(3e6);
end
% varsize declarations
if ~isempty(coder.target)
    if coder.internal.isConst(Ncfs) || ~coder.internal.eml_option_eq('UseMalloc','VariableSizeArrays')
        nnallvalues = coder.nullcopy(zeros(1,Ncfs,'like',x));
        prevnnallvalues = coder.nullcopy(zeros(1,UB,'like',x));    
        coder.varsize('leaders',[1 UB],[0 1]);
        leaderData = coder.nullcopy(zeros(1,UB,'like',x));
        coder.varsize('values',[1 UB],[0 1]);
        coder.varsize('nnallvalues',[1 UB],[0 1]);
        coder.varsize('prevnnallvalues',[1 UB],[0 1]);
        coder.varsize('approxcoefs',[1 UB],[0 1]);
        coder.varsize('details',[1 UB],[0 1]);
        coder.varsize('Abscoefs',[1 UB],[0 1]);
        coder.varsize('matneighbors',[3 UB],[0 1]);
        coder.varsize('matprevval',[3 UB],[0 1]);        
    else        
        nnallvalues = coder.nullcopy(zeros(1,UB,'like',x));
        prevnnallvalues = coder.nullcopy(zeros(1,UB,'like',x));
        leaderData = coder.nullcopy(zeros(1,UB,'like',x));
        coder.varsize('values',[1 UB],[0 1]);
        coder.varsize('nnallvalues',[1 UB],[0 1]);
        coder.varsize('prevnnallvalues',[1 UB],[0 1]);
        coder.varsize('approxcoefs',[1 UB],[0 1]);
        coder.varsize('details',[1 UB],[0 1]);
        coder.varsize('Abscoefs',[1 UB],[0 1]);
        coder.varsize('matneighbors',[3 UB],[0 1]);
        coder.varsize('matprevval',[3 UB],[0 1]);        
        coder.varsize('matneighbors',[3 UB],[0 1]);
        coder.varsize('leaderData',[1 UB],[0 1]);        
    end
    leaders = repmat({leaderData},1,Nlevels);
end

ncount = zeros(Nlevels,1);
if ~isempty(coder.target)
    coder.varsize('xtmp',[1 UB],[0 1]);
end
xtmp = x;
% Begin transform
for jj = 1:Nlevels
    nj = numel(xtmp);
    approxcoefs = conv(xtmp,Lo);
    assert(length(approxcoefs) == length(xtmp)+nwav-1);
    details = conv(xtmp,Hi);
    assert(length(details) == length(xtmp)+nwav-1);
    Na = length(xtmp)+nwav-1;
    % Length of detail coefficients matches approximation coefficients
    Nd = Na;
    % Set any NaNs equal to Infs
    approxcoefs(isnan(approxcoefs)) = Inf(1,1,'like',xtmp);
    details(isnan(details)) = Inf(1,1,'like',xtmp);
    % Remove boundary coefficients from scaling and wavelet coefficients
    if ~isempty(coder.target) && ~coder.internal.eml_option_eq('UseMalloc','VariableSizeArrays')
        % boundary coefficients at beginning 
        for ii = 1:nvalid-1
            approxcoefs(ii) = Inf(1,1,'like',xtmp);
            details(ii) = Inf(1,1,'like',xtmp);
        end
        % boundary coefficients at end 
        for ii = nj+1:Na
            approxcoefs(ii) = Inf(1,1,'like',xtmp);
            details(ii) = Inf(1,1,'like',xtmp);
        end
    else
        approxcoefs([1:nvalid-1 nj+1:Na]) = Inf(1,1,'like',xtmp);
        details([1:nvalid-1 nj+1:Nd]) = Inf(1,1,'like',xtmp);
    end
    % Decimate
    approxcoefs = approxcoefs((1:2:nj)+nwav-1);
    details = details((1:2:nj)+x0-1);
    % Replace data with approximation coefficients
    xtmp = approxcoefs;       
    % Use L1 normalization of coefficients
    Abscoefs = abs(details);
    Abscoefs = Abscoefs*2^(-jj/2);  
    Nd = length(Abscoefs);   
    % Determine wavelet leaders
    if jj == 1        
        %compute and store leaders
        nnallvalues = Abscoefs; %#ok<*AGROW>
        matneighbors = [Abscoefs(1,1:Nd-2) ; Abscoefs(1,2:Nd-1); Abscoefs(1,3:Nd)];
        neighbors = max(matneighbors);
        idxfinite = isfinite(neighbors);
        % Determine leaders
        leaders{jj} = neighbors(idxfinite);
        % How many wavelet leaders do we have at a given level
        ncount(jj) = numel(leaders{jj});             
    else
        nc = floor(numel(nnallvalues)/2);
        assert(nc <= n);
        matprevvall = [Abscoefs(1,1:nc); prevnnallvalues(1,1:2:2*nc); prevnnallvalues(1,2:2:2*nc)];
        nnallvalues = max(matprevvall);
        matneighbors = [nnallvalues(1,1:nc-2) ; nnallvalues(1,2:nc-1); nnallvalues(1,3:nc)];
        neighbors = max(matneighbors);
        idxfinite = isfinite(neighbors);
        leaders{jj} = neighbors(idxfinite);
        ncount(jj) = numel(leaders{jj});               
    end
    prevnnallvalues = nnallvalues;
end
coder.internal.errorIf(any(ncount == 0),'Wavelet:mfa:EmptyCount');
scales = 2.^(1:Nlevels);








