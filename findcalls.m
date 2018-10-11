function findcalls(wd,fs,varargin)
%% Find, extract, and save social vocalizations from audio recordings
% INPUT:
% wd: working directory, file path which contains files to be analyzed
% fs: sampling rate of recorded files
% Optional: 
% fileType: 'wav' or 'mat' to indicate format of recordings
% anal_dir: file path to save extracted files
%
% WARNING: Many individual parameters to be tuned!
%
% OUTPUT:
% Individual calls will be saved off in individual file in anal_dir named 
% with the original files's name and appeneded with a string '_Call_XXX'
% where XXX is the three digit number of calls within that file.
%%

pnames = {'fileType', 'outputDir', 'debug','params','dataVarName','high_filter_cutoff','low_cenv_filter_cutoff','adaptiveThreshold','thresh'};
dflts  = {'wav', fullfile(wd, 'Analyzed_auto'),false ,[],'recsGroup',2000,1000,false,0.75e-3};
[fileType,outputDir,debugFlag,params,dataVarName,high_filter_cutoff,low_cenv_filter_cutoff,adaptiveThreshold,thresh] = internal.stats.parseArgs(pnames,dflts,varargin{:});

if adaptiveThreshold
    n_std_thresh = 10;
    thresh = [];
end

if isempty(params)
    % (t>=calllength && H>rmsthresh && powerRatio<powerRatioThresh && wEnt<wEntThresh)
    callLength = 0.01; % in s
    mergethresh = 5e-3; % in s
    rmsthresh=0;
    powerRatioThresh = Inf;
    wEntThresh = Inf;
    
    params = struct('thresh',thresh,'callLength',callLength,'rmsthresh',rmsthresh,'mergethresh',mergethresh,'powerRatioThresh',powerRatioThresh,'wEntThresh',wEntThresh,'fs',fs);
end

rec_files = dir(fullfile(wd,['*.' fileType]));

if ~isfolder(outputDir)
    mkdir(outputDir);
end
n_files = length(rec_files);

for fln = 1:n_files
    filename=rec_files(fln).name;
    disp(filename)
    disp(['Analyzing file: ' num2str(fln) ' of ' num2str(n_files)])
    disp('...')
    
    switch fileType
        case 'wav'
            data_raw = audioread(fullfile(wd, filename));
        case 'mat'
            data_raw = load(fullfile(wd, filename));
            if isfield(data_raw,'analyzed')
               if data_raw.analyzed
                   continue
               end
            end
            data_raw = data_raw.(dataVarName);
    end
    
    [high_b, high_a] = butter(2,2*high_filter_cutoff/fs,'high');
    data = filtfilt(high_b,high_a,data_raw);
    if adaptiveThreshold
        data_round10 = 10*floor(length(data)/10);
        reshape_data_MA = reshape(data(1:data_round10),[],10);
        thresh = n_std_thresh*min(std(reshape_data_MA));
    end
    hilbenv = abs(hilbert(data));			
    [b,a] = butter(2,2*low_cenv_filter_cutoff/fs);			
    senv = filtfilt(b,a,hilbenv);
    allWins = calcWins(senv,params(1));
    [isCall,wins] = findCall(allWins,data_raw,params(1));
    
    if length(params)==2 && any(isCall)
        allWins = calcWins(senv,params(2));
        [isCall,wins] = findCall(allWins,data_raw,params(2));
    end
    
    if debugFlag
        if any(isCall)
            cla
            hold on
            plot(data);
            plot(allWins',max(data)*ones(2,size(allWins,1)));
            sound(data_raw,min([fs,200e3]));
            
            for w = find(isCall)
                callpos = wins(w,:);
                cut = data_raw(callpos(1):callpos(2));
                plot(callpos',max(data)*ones(2,1),'LineWidth',5);
                sound(cut,min(fs,200e3));
                keyboard;
            end
            plot(get(gca,'xlim'),[thresh thresh],'k')
        end
    else
        file_callcount = 0;
        for w = find(isCall)
            callpos = wins(w,:);
            cut = data_raw(callpos(1):callpos(2));
            save(fullfile(outputDir, [filename(1:end-4) '_Call_' sprintf('%03d',file_callcount) '.mat']),'cut','callpos','fs');
            file_callcount = file_callcount + 1;
        end
        if strcmp(fileType,'mat')
            save(fullfile(wd, filename),'analyzed','-append');
        end
    end
    
    
end

end

function wins = calcWins(senv,params)

thresh = params.thresh;
fs = params.fs;
mergethresh = params.mergethresh;

if length(unique(sign(senv-thresh))) > 1
    thresh_up = find(diff(sign(senv-thresh))==2);
    thresh_down = find(diff(sign(senv-thresh))==-2);
    if isempty(thresh_up) || isempty(thresh_down) % check that threhold crossing happens within file and not at beginning or end
        wins = [];
    else
        if length(thresh_up)==length(thresh_down) && thresh_up(1)<thresh_down(1)
            wins = [thresh_up,thresh_down];
        elseif length(thresh_up)>length(thresh_down)
            wins = [thresh_up(1:end-1),thresh_down];
        elseif length(thresh_down)>length(thresh_up)
            wins = [thresh_up,thresh_down(2:end)];
        elseif length(thresh_up)==length(thresh_down) && thresh_up(1)>thresh_down(1)
            wins = [thresh_up(1:end-1),thresh_down(2:end)];
        end
        if ~isempty(wins)
            wins = merge_wins(wins,fs,mergethresh);
        end
    end
else
    wins = [];
end

end

function [isCall,wins] = findCall(wins,data,params)

if isempty(wins)
    isCall = [];
    wins = [];
    return
end

calllength = params.callLength;
rmsthresh = params.rmsthresh;
mergethresh = params.mergethresh;
powerRatioThresh = params.powerRatioThresh;
wEntThresh = params.wEntThresh;
fs = params.fs;

wins = merge_wins(wins,fs,mergethresh);
isCall = false(1,size(wins,1));


for w = 1:size(wins,1)
    callpos = wins(w,:);
    cut = data(callpos(1):callpos(2));
    t=(length(cut)/fs);
    H=rms(cut);
    wEnt = weinerEntropy(cut);
    powerRatio = bandpower(cut,fs,[0 2e3])/bandpower(cut,fs,[3e3 5e3]);
    if (t>=calllength && H>rmsthresh && powerRatio<powerRatioThresh && wEnt<wEntThresh)
        isCall(w) = true;
    end
end

end

function WE = weinerEntropy(sig)

L = size(sig,1);

nfft = 2^nextpow2(L);
AFFT = fft(sig,nfft)./L;
AFFT = 2*abs(AFFT(1:nfft/2+1,:));
WE = exp(sum(log(AFFT)) ./ length(AFFT)) ./ mean(AFFT);

end

function [AFFT, F] = afft(sig,fs,nfft)

if size(sig,1) == 1
    sig = sig';
end

L = size(sig,1);

if nargin < 3 || isempty(nfft)
    nfft = 2^nextpow2(L);
end

F = fs/2*linspace(0,1,nfft/2+1);
AFFT = fft(sig,nfft)./L;
AFFT = 2*abs(AFFT(1:nfft/2+1,:));

end
