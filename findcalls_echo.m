function findcalls_echo(wd,fs,varargin)
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
    anal_dir = [wd 'Analyzed_auto_echo' filesep];
elseif nargin == 3
    fileType = varargin{1};
    anal_dir = [wd 'Analyzed_auto_echo' filesep];
elseif nargin == 4
    fileType = varargin{1};
    anal_dir = varargin{2};
end


debug = false;
low_filter_cutoff = [1e3 2e3];
high_filter_cutoff = [40e3 50e3];
cenv_filter_cutoff = [1e3 2e3];
n_std_thresh_high = 10;
n_std_thresh_low = 10;
quiet_file_thresh = 1e-3;
echo_separation = 1e-3 * [10 40]; %in sec
durthresh_max= 1e-2;
echo_peak_offset = 1e-3*[3 3];
n_baseline_chunks = 10;
smooth_span = 4e-3*fs;
min_peak_prominence = 1e-5;
remove_isolated_calls = false;

[hpFilt(1,:), hpFilt(2,:)] = ellip(5,5,60,2*high_filter_cutoff(1)/fs,'high');
[lpFilt(1,:), lpFilt(2,:)] = ellip(5,5,60,2*low_filter_cutoff(2)/fs,'low');
[cenvFilt(1,:), cenvFilt(2,:)] = ellip(5,5,60,2*cenv_filter_cutoff(2)/fs,'low');

rec_files = dir([wd '*.' fileType]);

if ~isfolder(anal_dir)
    mkdir(anal_dir);
end

if debug
    h1 = figure;
    h2 = figure;
end
n_files = length(rec_files);
total_callcount=0;
for fln = 1:n_files
    file_callcount = 0;
    filename=rec_files(fln).name;
    disp(filename)
    disp(['Analyzing file: ' num2str(fln) ' of ' num2str(n_files)])
    disp('...')
    data_raw = audioread([wd filename]);
    if debug
        figure(h1)
        plot((1:length(data_raw))/fs,data_raw);
        drawnow
    end
    
    if ~any(data_raw > quiet_file_thresh)
       continue 
    end
    
    senv_high = get_senv(data_raw,hpFilt,cenvFilt);
    senv_low = get_senv(data_raw,lpFilt,[]);
    thresh_high = get_thresh(senv_high,n_baseline_chunks,n_std_thresh_high);
    thresh_low = get_thresh(senv_low,n_baseline_chunks,n_std_thresh_low);
    [pks, locs] = findpeaks(smooth(senv_high,smooth_span),fs,'MinPeakHeight',thresh_high,'MaxPeakWidth',durthresh_max,'MinPeakDistance',echo_separation(1),'MinPeakProminence',min_peak_prominence);
    echo_peaks = false(1,length(pks));
    checked_peaks = false(1,length(pks));
    low_thresh_check = false(1,length(pks));
    echo_length_check = false(1,length(pks));
    iscall = false(1,length(pks));
    wins = zeros(length(pks),2);
    while(~all(checked_peaks))
        w = find(~checked_peaks,1);
        echo_peak = find(locs > locs(w) + echo_separation(1) &  locs < locs(w) + echo_separation(2), 1);
        if length(echo_peak) == 1
            call_win = round([max([(locs(w)-echo_peak_offset(1)) 1/fs]) min([(locs(echo_peak)+echo_peak_offset(2)) length(data_raw)/fs])]*fs);
            if all(senv_low(call_win(1):call_win(2))<thresh_low) 
                echo_peaks([w, echo_peak]) = true;
                checked_peaks([w, echo_peak]) = true;
                wins(w,:) = call_win;
                iscall(w) = true;
                file_callcount = file_callcount + 1;
                total_callcount = total_callcount + 1;
                disp(['Call count: ' num2str(total_callcount)])
            else
                low_thresh_check([w, echo_peak]) = true;
                checked_peaks([w, echo_peak]) = true;
            end
        else
            echo_length_check([w, echo_peak]) = true;
            checked_peaks([w, echo_peak]) = true;
        end
    end
    all_calls = wins;
    wins = wins(wins(:,1)~=0,:);
    if remove_isolated_calls
        if file_callcount >= 2
            call_sep = [inf; wins(2:end,1) - wins(1:end-1,2); inf];
            isolated_calls = zeros(1,file_callcount);
            for call = 1:file_callcount
                isolated_calls(call) = ~any(call_sep(call:call+1) <  min_call_sep);
            end
            wins = wins(~isolated_calls,:);
            if sum(~isolated_calls) > 2
                file_callcount = file_callcount - sum(isolated_calls);
                total_callcount = total_callcount - sum(isolated_calls);
            else
                total_callcount = total_callcount - file_callcount;
                file_callcount = 0;
            end
            
        else
            total_callcount = total_callcount - file_callcount;
            file_callcount = 0;
        end
    end
    wins = merge_wins(wins,fs,20e-3);
    cutcalls = cell(1,size(wins,1));
    for w = 1:size(wins,1)
       cutcalls{w} = data_raw(wins(w,1):wins(w,2));
    end
    if debug && file_callcount
        for w = 1:size(wins,1)
            figure(h1);
            clf
            subplot(1,2,1);
            sound(cutcalls{w},10e3);
            spectrogram(cutcalls{w},kaiser(512,6),510,2048,fs,'yaxis');
            subplot(1,2,2);
            hold on
            plot(senv_high(wins(w,1):wins(w,2)));
            plot([0 diff(wins(w,:))],[thresh_high thresh_high],'g')
            plot(senv_low(wins(w,1):wins(w,2)),'r');
            plot(fs*[0 echo_peak_offset(1)+durthresh_max],[thresh_low thresh_low],'k')
            plot(diff(wins(w,:)) - fs*[echo_peak_offset(2)+durthresh_max 0],[thresh_low thresh_low],'k')
            keyboard;
        end
        figure(h2);
        hold on
        t = 1e3*(1:length(data_raw))/fs;
        plot(t,data_raw)
        plot((1e3*wins/fs)',0.1*ones(2,size(wins,1)),'LineWidth',12);
        plot(t,senv_high,'r')
        plot([0 t(end)],[thresh_high thresh_high],'g')
        if any(echo_peaks)
            plot(1e3*(locs(echo_peaks)),pks(echo_peaks).*ones(sum(echo_peaks),1),'ro')
        end
        keyboard;
        cla;
    elseif file_callcount
        for w = 1:size(wins,1)
            cut = cutcalls{w};
            callpos = wins(w,:);
            save([anal_dir filename(1:end-4) '_Echo_' sprintf('%03d',w-1) '.mat'],'cut','callpos','fs');
        end
    end
end
end

function baseline = get_thresh(data,n_chunks,n_std_thresh)

data_length_round = n_chunks*floor(length(data)/n_chunks);
data_reshape = reshape(data(1:data_length_round),[],n_chunks);
[min_std, min_std_idx] = min(std(data_reshape));
baseline = n_std_thresh*min_std + mean(data_reshape(:,min_std_idx));

end

function senv = get_senv(data,dataFilt,cenvFilter)

data_filt = filtfilt(dataFilt(1,:),dataFilt(2,:),data);
hilbenv = abs(hilbert(data_filt));
if ~isempty(cenvFilter)
    senv = filtfilt(cenvFilter(1,:),cenvFilter(2,:),hilbenv);
else
    senv = hilbenv;
end

end
