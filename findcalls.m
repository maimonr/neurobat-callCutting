function [fcns,call_cut_params] = findcalls(data_raw,wd,fs,varargin) 

pnames = {'fileType', 'outputDir', 'debug','params','dataVarName','envelope_win_size','adaptiveThreshold','thresh','n_std_thresh','filename','next_filename','plot_offset_s'};
dflts  = {'wav', fullfile(wd, 'Analyzed_auto'),false ,[],'recsGroup',1e-3,false,0.75e-3,10,'cut_call',[],3};
[fileType,outputDir,debugFlag,call_cut_params,dataVarName,envelope_win_size,adaptiveThreshold,thresh,n_std_thresh,filename,next_filename,plot_offset_s] = internal.stats.parseArgs(pnames,dflts,varargin{:});

allParams = struct('fileType',fileType,'outputDir',outputDir,'debugFlag',debugFlag,...
    'dataVarName',dataVarName,'adaptiveThreshold',adaptiveThreshold,'n_std_thresh',n_std_thresh);

if adaptiveThreshold
    thresh = [];
end

if isempty(call_cut_params)
    % (t>=calllength && H>rmsthresh && powerRatio<powerRatioThresh && wEnt<wEntThresh)
    callLength = 0.015; % in s
    mergethresh = 5e-3; % in s
    rmsthresh=0;
    powerRatioThresh = Inf;
    wEntThresh = Inf;
    
    call_cut_params = struct('thresh',thresh,'callLength',callLength,...
        'rmsthresh',rmsthresh,'mergethresh',mergethresh,...
        'powerRatioThresh',powerRatioThresh,'wEntThresh',wEntThresh,'fs',fs,...
        'n_std_thresh',n_std_thresh,'envelope_win_size',envelope_win_size);
end

if nargout > 0
    fcns = localfunctions;
    return
end

senv = get_envelope(data_raw,call_cut_params);

if adaptiveThreshold % if requested set threshold to the minimum standard deviation across 10 equal sized chunks across file
    data_round10 = 10*floor(length(data_raw)/10);
    reshape_data_MA = reshape(data(1:data_round10),[],10);
    thresh = n_std_thresh*min(std(reshape_data_MA));
    call_cut_params.thresh = thresh;
end

[allWins,multi_file_wins] = calcWins(senv,call_cut_params(1)); % generate potential windows containing calls
multi_file_wins = multi_file_wins{1};

if ~isempty(multi_file_wins) && ~isempty(next_filename) % check for calls overlapping between 2 .wav files
    % get succeeding file's data
    next_data_raw = load_audio_data(wd,next_filename,fileType,dataVarName);
    next_senv = get_envelope(next_data_raw,call_cut_params);
    next_senv = [senv(multi_file_wins(1):end); next_senv];
    % get end sample of call window in next file
    [next_file_wins, next_file_multi_wins] = calcWins(next_senv,call_cut_params(1));
    if ~isempty(next_file_multi_wins{2})
        next_file_wins = vertcat(next_file_multi_wins{2},next_file_wins);
    end
    next_file_wins = merge_wins(next_file_wins,fs,mergethresh);
    next_file_wins = next_file_wins(1,:) - diff(multi_file_wins)+1;
    % concatenate data across files
    multi_file_data_raw = [data_raw(multi_file_wins(1):end); next_data_raw(1:next_file_wins(2))];
    multi_file_win = [1 length(multi_file_data_raw)];
    isCall = findCall(multi_file_win,multi_file_data_raw,call_cut_params(1));
    
    if isCall && ~debugFlag
        callpos = [multi_file_wins(1) next_file_wins(1,2)];
        cut = multi_file_data_raw;
        save(fullfile(outputDir, [filename '_' next_filename(1:end-4) '_Call_' sprintf('%03d',1) '.mat']),'cut','callpos','fs','allParams','call_cut_params');
    end
end

[isCall,wins] = findCall(allWins,data_raw,call_cut_params(1)); % merge windows and assess if they contain calls

if length(call_cut_params)==2 && any(isCall) % do a second pass if requested
    allWins = calcWins(senv,call_cut_params(2));
    [isCall,wins] = findCall(allWins,data_raw,call_cut_params(2));
end

if debugFlag
    if any(isCall)
        hold on
        for w = find(isCall)
            %%
            cla
            clear sound
            callpos = wins(w,:);
            call_pos_win = callpos + (fs*[-plot_offset_s plot_offset_s]);
            call_pos_win(1) = max(1,call_pos_win(1)); call_pos_win(2) = min(length(data_raw),call_pos_win(2));
            call_pos_idx = call_pos_win(1):call_pos_win(2);
            
            plot(allWins'/fs,max(data_raw)*ones(2,size(allWins,1)));
            plot(callpos'/fs,0.9*max(data_raw(call_pos_idx))*ones(2,1),'LineWidth',5);
            plot(call_pos_idx/fs,data_raw(call_pos_idx));
            xlim(call_pos_win/fs);
            ylim([min(data_raw(call_pos_idx)) max(data_raw(call_pos_idx))])
            sound(data_raw(call_pos_idx),min(fs,200e3));
            %%
            keyboard;
        end
        plot(get(gca,'xlim'),[thresh thresh],'k')
    end
else
    file_callcount = 0;
    for w = find(isCall) % save off calls and detection parameters
        callpos = wins(w,:);
        cut = data_raw(callpos(1):callpos(2));
        save(fullfile(outputDir, [filename '_Call_' sprintf('%03d',file_callcount) '.mat']),'cut','callpos','fs','allParams','call_cut_params');
        file_callcount = file_callcount + 1;
    end
    if strcmp(fileType,'mat')
        save(fullfile(wd, filename),'analyzed','-append');
    end
end


end

function senv = get_envelope(data_raw,call_cut_params)

windowLength = round(call_cut_params(1).fs*call_cut_params(1).envelope_win_size);
senv = envelope(data_raw,windowLength,'rms');

end

function data_raw = load_audio_data(wd,filename,fileType,dataVarName)

switch fileType
    case 'wav'
        data_raw = audioread(fullfile(wd, filename));
    case 'mat'
        data_raw = load(fullfile(wd, filename));
        if isfield(data_raw,'analyzed')
            if data_raw.analyzed
                data_raw = [];
                return
            end
        end
        data_raw = data_raw.(dataVarName);
end

end

function [wins,multi_file_wins] = calcWins(senv,params)

thresh = params.thresh;
fs = params.fs;
mergethresh = params.mergethresh;
multi_file_wins = cell(1,2);

if length(unique(sign(senv-thresh))) > 1
    thresh_up = find(diff(sign(senv-thresh))==2);
    thresh_down = find(diff(sign(senv-thresh))==-2);
    if isempty(thresh_up) && isempty(thresh_down) 
        wins = [];
    else
        if length(thresh_up)==length(thresh_down) && thresh_up(1)<thresh_down(1) % (A) all windows contained in current .wav file
            wins = [thresh_up,thresh_down];
        elseif length(thresh_up)>length(thresh_down) % (B) one window starts in this file and continues until the end of the file
            wins = [thresh_up(1:end-1),thresh_down];
            multi_file_wins{1} = [thresh_up(end) length(senv)];
        elseif length(thresh_down)>length(thresh_up) % (C) one window starts at the first sample of this file
            wins = [thresh_up,thresh_down(2:end)];
            multi_file_wins{2} = [1 thresh_down(1)]; 
        elseif length(thresh_up)==length(thresh_down) && thresh_up(1)>thresh_down(1) % (B) AND (C)
            wins = [thresh_up(1:end-1),thresh_down(2:end)];
            multi_file_wins{1} = [thresh_up(end) length(senv)];
            multi_file_wins{2} = [1 thresh_down(1)];
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
