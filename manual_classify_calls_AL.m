function manual_classify_calls_AL(expDate,call_str)
al_fs = 50e3;
addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
dataDir = 'D:\acoustic_recording\';
T = readtable([dataDir 'recording_logs.csv']);
T = T(T.Date==expDate,:);
logger_base_dir = 'Z:\users\Maimon\acoustic_recording\video_and_AL\';
audio_base_dir = 'Z:\users\Maimon\acoustic_recording\audio\';

date_str_format = 'mmddyyyy';
dateStr = datestr(expDate,date_str_format);

audio_dir = [audio_base_dir dateStr filesep 'audio\ch1\'];

switch call_str
    case 'Call'
        callDir = [audio_dir 'Analyzed_auto' filesep];
    case 'Echo'
        callDir = [audio_dir 'Analyzed_auto_echo' filesep];
end

audio2nlg = load([audio_dir 'audio2nlg_fit.mat']);

logger_nums = T{1,4:2:14};
logger_nums = logger_nums(~isnan(logger_nums));
logger_nums = setdiff(logger_nums,str2double(T.malfunction_loggers{1}));
loggerDir = [logger_base_dir dateStr '\audiologgers\'];
for logger_k = 1:length(logger_nums)
    data_fname = dir([loggerDir 'logger' num2str(logger_nums(logger_k)) filesep 'extracted_data' filesep '*CSC0.mat']);
    tsData(logger_k) = load(fullfile(data_fname.folder,data_fname.name));
end

specParams = struct('colormap_name','hot','spec_win_size',512,'spec_overlap_size',500,'spec_nfft',1024,'fs',al_fs,'spec_ylims',[0 10],'spec_caxis_factor',0.1);
orig_rec_plot_win = 5;
recVar = 'cut';
recVar_full = 'recsGroup';

wav_mat_file = input('wav (1) or mat (2) file?');

callFiles = dir([callDir '*' call_str '*.mat']);
callNums = 1:length(callFiles);

if exist([callDir 'current_classify_file_number.mat'],'file')
    f = load([callDir 'current_classify_file_number.mat']);
    fNum = find(callNums == abs(min(callNums-f.fNum)));
else
    fNum = 1;
end
cum_samples = [0 cumsum(audio2nlg.total_samples_by_file)];
nCalls = length(callFiles);
nLogger = length(tsData);

for c = fNum:nCalls
    
    s = load([callDir callFiles(c).name]);
    data = s.(recVar);
    fs = min(s.fs,200e3);
    sound(data,fs);
    origRec_fName_split = strsplit(callFiles(c).name,'_');
    wav_f_num = str2double(origRec_fName_split{2});
    
    if wav_mat_file == 1
        if length(origRec_fName_split) > 4
            origRec_fName = fullfile(audio_dir, [strjoin(origRec_fName_split(1:2),'_') '.WAV']);
            sample_offset = audio2nlg.total_samples_by_file(wav_f_num) - s.callpos(1);
            data_full_append = s.cut(sample_offset:end);
            dataFull = vertcat(audioread(origRec_fName),data_full_append);
            s.callpos = [s.callpos(1) length(dataFull)];
        else
            origRec_fName = fullfile(audio_dir, [strjoin(origRec_fName_split(1:end-2),'_') '.WAV']);
            dataFull = audioread(origRec_fName);
        end
    elseif wav_mat_file == 2
        origRec_fName = fullfile(audio_dir, [strjoin(origRec_fName_split(1:end-2),'_') '.mat']);
        d = load(origRec_fName);
        dataFull = d.(recVar_full);
    end
    
    callNum = strsplit(callFiles(c).name,{'_','.'});
    callNum = str2double(callNum{end-1});
    %%
    subplot(nLogger+1,1,nLogger+1)
    hold on
    if callNum == 0
        cla
        plot((1:length(dataFull))/s.fs,dataFull,'k');
    end
    plot((s.callpos(end,1)+(0:length(data)-1))/s.fs,data);
    
    if length(dataFull)/s.fs > orig_rec_plot_win
       xlim([s.callpos(end,1)/s.fs - orig_rec_plot_win/2 s.callpos(end,1)/s.fs + orig_rec_plot_win/2])
    end
    
    corrected_callpos = (1e3*(s.callpos + cum_samples(wav_f_num))/s.fs);
    corrected_callpos = avi2nlg_time(audio2nlg,corrected_callpos);
    
    call_ts_nlg = (corrected_callpos + audio2nlg.first_nlg_pulse_time)*1e3;
    
    piezo_chunk = cell(1,nLogger);
    specParams = struct('colormap_name','hot','spec_win_size',512,'spec_overlap_size',500,'spec_nfft',1024,'fs',al_fs,'spec_ylims',[0 10],'spec_caxis_factor',0.1);
    logger_k = 1;
    call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
        tsData(logger_k).Indices_of_first_and_last_samples(:,1)',...
        tsData(logger_k).Timestamps_of_first_samples_usec,...
        tsData(logger_k).Sampling_period_usec_Logger);
    if isempty(call_ts_nlg_sample)
        break
    end
    
    for logger_k = 1:nLogger
        subplot(nLogger+1,1,logger_k)
        cla
        call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
            tsData(logger_k).Indices_of_first_and_last_samples(:,1)',...
            tsData(logger_k).Timestamps_of_first_samples_usec,...
            tsData(logger_k).Sampling_period_usec_Logger);
        if isempty(call_ts_nlg_sample)
            break
        end
        piezo_chunk_idx = call_ts_nlg_sample(1):call_ts_nlg_sample(end);
        piezo_chunk{logger_k} = double(tsData(logger_k).AD_count_int16(piezo_chunk_idx));
        
        plot_call_spectrogram(piezo_chunk{logger_k},gca,specParams);
    end
    %%
    display(callFiles(c).name);
    
    repeat = 1;
    repeat_k = 1;
    while repeat
        class = input('call?','s');
        switch class
            case '1'
                repeat = 0;
                s.noise = false;
                save([callDir callFiles(c).name],'-struct','s');
            case '0'
                repeat = 0;
                s.noise = true;
                save([callDir callFiles(c).name],'-struct','s');
            case 'stop'
                fNum = callNums(c);
                save([callDir 'current_classify_file_number.mat'],'fNum');
            case 'pause'
                keyboard
            otherwise
                pause(0.1);
                if repeat_k < 3
                    sound(data,fs/(2*repeat_k));
                else
                    startIdx = max(1,s.callpos(end,1) - (orig_rec_plot_win/2)*s.fs);
                    endIdx = min(length(dataFull),s.callpos(end,1) + (orig_rec_plot_win/2)*s.fs);
                    sound(dataFull(startIdx:endIdx),fs);
                    repeat_k = 1;
                end
                repeat_k = repeat_k + 1;
        end
    end
end

save_call_data = input('build and save cut call data file?');
if save_call_data
    cut_call_data = get_corrected_call_times(audio_dir,callDir,call_str);
    
    switch call_str
        case 'Call'
            save([audio_dir 'cut_call_data.mat'],'cut_call_data');
            
        case 'Echo'
            save([audio_dir 'cut_echo_data.mat'],'cut_call_data');
            
    end
end

end