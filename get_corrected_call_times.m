function cut_call_data = get_corrected_call_times(audioDir,analyzed_audio_dir,call_str)

audio2nlg = load(fullfile(audioDir, 'audio2nlg_fit.mat')); % load correction between avisoft and NLG time data

wav_files = dir(fullfile(audioDir, '*.wav')); % load raw audio recordings
% sort wav files according to numbering
wav_files_name = {wav_files.name};
wav_file_nums = cellfun(@(x) str2double(regexp(x,'(?<=_)\d+(?=.WAV)','match','ignorecase')), wav_files_name); 
[~, sort_wav_files_idx] = sort(wav_file_nums);

% load cut call files
cut_call_files = dir([analyzed_audio_dir 'T*' call_str '*.mat']);
n_cut_call_files = length(cut_call_files);

% get experiment date as datetime
exp_day_str = 'neurologger_recording';
dateFormat = 'yyyyMMdd';
idx = strfind(audioDir,exp_day_str) + length(exp_day_str);
expDay = datetime(audioDir(idx:idx+length(dateFormat)-1),'InputFormat',dateFormat);

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
start_f_num = wav_file_nums(sort_wav_files_idx(1));

% calculate cumulative samples from beginning of recording to place call in
% time relative to entire recording
cum_samples = [0 cumsum(audio2nlg.total_samples_by_file)];
if any(isnan(cum_samples))
    goodCorrection = 0;
else
    goodCorrection = 1;
end
fs_wav = 250e3; % add in minor sampling rate correction

for call_f = 1:n_cut_call_files
    s = load([analyzed_audio_dir cut_call_files(call_f).name]);
    % transfer all data from cut call file to the data structure we're
    % building here
    cut_call_fields = fieldnames(s); 
    for f = 1:length(cut_call_fields)
        cut_call_data(call_f).(cut_call_fields{f}) = s.(cut_call_fields{f});
    end
    f_name_split = strsplit(cut_call_files(call_f).name,{'_','.'});
    cut_call_data(call_f).fName = fullfile(audioDir, [strjoin(f_name_split(1:2),'_') '.WAV']);
    cut_call_data(call_f).f_num = str2double(f_name_split{2}) - start_f_num + 1;
    cut_call_data(call_f).expDay = expDay;
    cut_call_data(call_f).batNum = batNum;
    % calculate time in whole recordings, corrected from AVI/NLG clock
    % drift
    if goodCorrection
        cut_call_data(call_f).corrected_callpos = (1e3*(cut_call_data(call_f).callpos + cum_samples(cut_call_data(call_f).f_num))/fs_wav);
        cut_call_data(call_f).corrected_callpos = avi2nlg_time(audio2nlg,cut_call_data(call_f).corrected_callpos);
    else
        cut_call_data(call_f).corrected_callpos = nan(1,2);
    end
end

cut_call_fields = fieldnames(cut_call_data);
cut_call_fields = cut_call_fields(~strcmp(cut_call_fields,'noise'));
for f = 1:length(cut_call_fields)
    emptyFields = arrayfun(@(x) isempty(x.(cut_call_fields{f})), cut_call_data);
    assert(~any(emptyFields));
end

end