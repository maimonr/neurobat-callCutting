function findcalls_v5_mcr(wd,fs,varargin)

if nargin == 2
    wav_mat_file = input('wav (1) or mat (2) file?');
    if wav_mat_file == 1
        fileType = 'wav';
    elseif wav_mat_file == 2
        fileType = 'mat';
    else
        disp('Invalid input');
        return
    end
    anal_dir = [wd 'Analyzed_auto' filesep];
elseif nargin == 3
    fileType = varargin{1};
    anal_dir = [wd 'Analyzed_auto' filesep];
elseif nargin == 4
    fileType = varargin{1};
    anal_dir = varargin{2};
end

reanalyze = false;
dataVarName = 'recsGroup';
debugFlag = false;
callLength = 0.05; % in s
mergethresh = 5e-3; % in s
rmsthresh=0.015;
powerRatioThresh = 3;
wEntThresh = 0.25;
analyzed = true;

high_filter_cutoff = 2000;
low_cenv_filter_cutoff = 1000;
adaptiveThreshold = false;

if adaptiveThreshold
    n_std_thresh = 10;
else
    thresh = 1.2e-3;
end

params(1) = struct('thresh',thresh,'callLength',callLength,'rmsthresh',rmsthresh,'mergethresh',mergethresh,'powerRatioThresh',powerRatioThresh,'wEntThresh',wEntThresh,'fs',fs);
params(2) = params(1);
params(2).thresh = 0.75e-3;
params(2).callLength = 0.02;
params(2).rmsthresh = 0.001;
params(2).powerRatioThresh = 10;
params(2).wEntThresh = 0.5;

rec_files = dir([wd '*.' fileType]);

if ~isdir(anal_dir)
    mkdir(anal_dir);
end
n_files = length(rec_files);

for fln = 1:n_files
    filename=rec_files(fln).name;
    disp(filename)
    disp(['Analyzing file: ' num2str(fln) ' of ' num2str(n_files)])
    disp('...')
    
    switch fileType
        case 'wav'
            data_raw = audioread([wd filename]);
        case 'mat'
            data_raw = load([wd filename]);
            if isfield(data_raw,'analyzed') && ~reanalyze
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
    [wins,indx] = calcWins(senv,data,params(1));
    [isCall,wins] = findCall(wins,indx,data_raw,params(1));
    
    if any(isCall)
        [allWins,indx] = calcWins(senv,data,params(2));
        [isCall,wins] = findCall(allWins,indx,data_raw,params(2));
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
            save([anal_dir filename(1:end-4) '_Call_' sprintf('%03d',file_callcount) '.mat'],'cut','callpos','fs');
            file_callcount = file_callcount + 1;
        end
        if strcmp(fileType,'mat')
            save([wd filename],'analyzed','-append');
        end
    end
    
    
end

end

function [wins,indx] = calcWins(senv,data,params)

thresh = params.thresh;
callLength = params.callLength;
rmsthresh = params.rmsthresh;
fs = params.fs;
mergethresh = params.mergethresh;

if length(unique(sign(senv-thresh))) > 1
    thresh_up = find(diff(sign(senv-thresh))==2);
    thresh_down = find(diff(sign(senv-thresh))==-2);
    if isempty(thresh_up) || isempty(thresh_down) % check that threhold crossing happens within file and not at beginning or end
        wins = [];
        indx = [];
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
            win_cell = mat2cell(wins,ones(1,size(wins,1)),2);
            diffs = cellfun(@diff,win_cell)/fs;
            rms_win = cellfun(@(x) rms(data(x(1):x(2))),win_cell);
            indx=find((diffs>callLength) & (rms_win>rmsthresh));
        else
            indx = [];
        end
    end
else
    wins = [];
    indx = [];
end

end

function [isCall,wins] = findCall(wins,indx,data,params)

calllength = params.callLength;
rmsthresh = params.rmsthresh;
mergethresh = params.mergethresh;
powerRatioThresh = params.powerRatioThresh;
wEntThresh = params.wEntThresh;
fs = params.fs;

wins = merge_wins(wins(indx,:),fs,mergethresh);
isCall = false(1,size(wins,1));


for w = 1:size(wins,1)
    callpos = wins(w,:);
    cut = data(callpos(1):callpos(2));
    t=(length(cut)/fs);
    H=rms(cut);
    wEnt = weinerEntropy(cut);
    powerRatio = bandpower(cut,fs,[0 5e3])/bandpower(cut,fs,[5e3 10e3]);
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