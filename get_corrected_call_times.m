function cut_call_data = get_corrected_call_times(audioDir,analyzed_audio_dir,call_str,varargin)

pnames = {'audio_file_filter', 'recType','fs'};
dflts  = {'*', 'avisoft',250e3};
[audio_file_filter,recType,fs] = internal.stats.parseArgs(pnames,dflts,varargin{:});

audio2nlg = load(fullfile(audioDir, 'audio2nlg_fit.mat')); % load correction between avisoft and NLG time data

wav_files = dir(fullfile(audioDir, [audio_file_filter '.wav'])); % load raw audio recordings
% sort wav files according to numbering
wav_files_name = {wav_files.name};
wav_file_nums = cellfun(@(x) str2double(regexp(x,'(?<=_)\d+(?=.WAV)','match','ignorecase')), wav_files_name);
wav_file_nums = sort(wav_file_nums);

% load cut call files
cut_call_files = dir(fullfile(analyzed_audio_dir,['*_' call_str '_*.mat']));
n_cut_call_files = length(cut_call_files);

% get experiment date as datetime
exp_day_reg_exp_str = '\d{8}';
if contains(audioDir,'adult')
    dateFormat = 'MMddyyyy';
else
    dateFormat = 'yyyyMMdd';
end
exp_day_str = regexp(audioDir,exp_day_reg_exp_str,'match');
exp_day_str = exp_day_str{1};
expDay = datetime(exp_day_str,'InputFormat',dateFormat);

bat_num_str = 'bat';
bat_num_length = 5;
idx = strfind(audioDir,bat_num_str) + length(bat_num_str);
batNum = audioDir(idx:idx+bat_num_length-1);
if isempty(batNum)
    batNum = NaN;
end

% initialize data structure
cut_call_data = struct('callpos',{},'cut',{},'corrected_callpos',{},'f_num',{},'fName',{},'fs',{},'noise',{},'expDay',{},'batNum',{});
% determine which number raw .wav files started, in case not 1
start_f_num = wav_file_nums(1);

% calculate cumulative samples from beginning of recording to place call in
% time relative to entire recording
cum_samples = [0 cumsum(audio2nlg.total_samples_by_file)];
if any(isnan(cum_samples))
    goodCorrection = 0;
else
    goodCorrection = 1;
end

for call_f = 1:n_cut_call_files
    s = load(fullfile(analyzed_audio_dir, cut_call_files(call_f).name));
    % transfer all data from cut call file to the data structure we're
    % building here
    
    cut_call_fields = fieldnames(s);
    for f = 1:length(cut_call_fields)
        cut_call_data(call_f).(cut_call_fields{f}) = s.(cut_call_fields{f});
    end
    f_name_split = strsplit(cut_call_files(call_f).name,{'_','.'});
    
    
    switch recType
        
        case 'avisoft'
            
            if length(f_name_split) > 5
                fName = {fullfile(audioDir, [strjoin(f_name_split(1:2),'_') '.WAV']),...
                    fullfile(audioDir, [strjoin(f_name_split(3:4),'_') '.WAV'])};
                fNum = [str2double(f_name_split{2}) str2double(f_name_split{4})] - start_f_num + 1;
                if s.callpos(2) > audio2nlg.total_samples_by_file(fNum(1))
                    cut_call_data(call_f).callpos(2) = s.callpos(2) -  audio2nlg.total_samples_by_file(fNum(1));
                end
            else
                fName = fullfile(audioDir, [strjoin(f_name_split(1:2),'_') '.WAV']);
                fNum = str2double(f_name_split{2}) - start_f_num + 1;
            end
            
        case 'vocOperant'
            
            fName = fullfile(audioDir, [strjoin(f_name_split(1:6),'_') '.WAV']);
            fNum = str2double(f_name_split{6}) - start_f_num + 1;
            
    end
    
    cut_call_data(call_f).fName = fName;
    cut_call_data(call_f).f_num = fNum;
    cut_call_data(call_f).expDay = expDay;
    cut_call_data(call_f).batNum = batNum;
    % calculate time in whole recordings, corrected from AVI/NLG clock
    % drift
    if goodCorrection
        cut_call_data(call_f).corrected_callpos = (1e3*(cut_call_data(call_f).callpos + cum_samples(cut_call_data(call_f).f_num))/fs);
        cut_call_data(call_f).corrected_callpos = avi2nlg_time(audio2nlg,cut_call_data(call_f).corrected_callpos);
    else
        cut_call_data(call_f).corrected_callpos = nan(1,2);
    end
end
not_manually_classified_idx = arrayfun(@(x) isempty(x.noise),cut_call_data);
[cut_call_data(not_manually_classified_idx).noise] = deal(NaN);
cut_call_fields = fieldnames(cut_call_data);
for f = 1:length(cut_call_fields)
    emptyFields = arrayfun(@(x) isempty(x.(cut_call_fields{f})), cut_call_data);
    assert(~any(emptyFields));
end

all_f_num = arrayfun(@(x) x.f_num(1),cut_call_data);
[~,f_num_sort_idx] = sort(all_f_num);
cut_call_data = cut_call_data(f_num_sort_idx);

end