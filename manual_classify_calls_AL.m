function manual_classify_calls_AL(audio_dir,tsData,audio2nlg,varargin)
al_fs = 50e3;
addpath('C:\Users\phyllo\Documents\GitHub\SoundAnalysisBats\')
addpath('C:\Users\phyllo\Documents\Maimon\acoustic_recording\scripts\')
addpath('C:\Users\phyllo\Documents\GitHub\LoggerDataProcessing\')
addpath('C:\Users\phyllo\Documents\GitHub\neurobat-hardware-alignment\')

pnames = {'recType', 'classify_call_nums'};
dflts  = {'avisoft', []};
[expType,classify_call_nums] = internal.stats.parseArgs(pnames,dflts,varargin{:});

callDir = fullfile(audio_dir,'Analyzed_auto');

specParams = struct('colormap_name','hot','spec_win_size',512,'spec_overlap_size',500,'spec_nfft',512,'fs',al_fs,'spec_ylims',[0 10],'spec_caxis_factor',0.1);
orig_rec_plot_win = 5;
recVar = 'cut';

AL_call_offset_length = 1;
AL_call_offset = 1e3*0.5*AL_call_offset_length*[-1 1];
ds_factor = 10;

avi_fs_correction = 0;

callFiles = dir([callDir filesep '*Call*.mat']);
callNums = 1:length(callFiles);
nCalls = length(callFiles);

if strcmp(expType,'operant')
    wav_file_num_regexp_str = '_\d{1,2}_Call';
    wav_file_nums = arrayfun(@(x) regexp(x.name,wav_file_num_regexp_str,'match'),callFiles,'un',0);
    wav_file_nums = cellfun(@(x) str2double(strrep(strrep(x{1},'_',''),'Call','')),wav_file_nums);
    [~,wav_file_sort_idx] = sort(wav_file_nums);
    callFiles = callFiles(wav_file_sort_idx);
end

if isempty(classify_call_nums)
    if exist(fullfile(callDir,'current_classify_file_number.mat'),'file')
        f = load(fullfile(callDir,'current_classify_file_number.mat'));
        fNum = f.fNum;
    else
        fNum = 1;
    end
    classify_call_nums = fNum:nCalls;
end

manual_bat_ID_fName = fullfile(audio_dir,'manual_al_classify_batNum.mat');

if exist(manual_bat_ID_fName,'file')
    f = load(manual_bat_ID_fName);
    manual_al_classify_batNum = f.manual_al_classify_batNum;
else
    manual_al_classify_batNum = cell(1,nCalls);
    save(manual_bat_ID_fName,'manual_al_classify_batNum')
end

cum_samples = [0 cumsum(audio2nlg.total_samples_by_file)];
nLogger = length(tsData);
subplot_handles = cell(1,nLogger+1);

for logger_k = 1:nLogger
    subplot_handles{logger_k} = subplot(nLogger+1,1,logger_k);
    hold(subplot_handles{logger_k},'on');
    subplot_handles{logger_k}.XAxis.Visible = 'off';
end

subplot_handles{end} = subplot(nLogger+1,1,nLogger+1);
hold(subplot_handles{end},'on');

last_rec_name = [];

for c = classify_call_nums
    
    call_file_fName = fullfile(callDir,callFiles(c).name);
    s = load(call_file_fName);
    data = s.(recVar);
    fs = min(s.fs,200e3);
    avi_fs = s.fs - avi_fs_correction;
    sound(data,fs);
   
    origRec_fName_split = strsplit(callFiles(c).name,'_');
    
     if contains(callFiles(c).name,'VocTrigger')
         origRec_fName = fullfile(audio_dir, [strjoin(origRec_fName_split(1:6),'_') '.WAV']);
         split_file_flag = false;
         if length(origRec_fName_split) > 10
             noiseStatus = true;
             manual_al_classify_batNum{c} = 'noise';
             save_bat_num(s,noiseStatus,manual_al_classify_batNum,manual_bat_ID_fName,call_file_fName)
             continue
         end
     else
         if length(origRec_fName_split) > 4
             split_file_flag = true;
             origRec_fName = fullfile(audio_dir, [strjoin(origRec_fName_split(1:2),'_') '.WAV']);
         else
             split_file_flag = false;
             origRec_fName = fullfile(audio_dir, [strjoin(origRec_fName_split(1:end-2),'_') '.WAV']);
         end
         
     end
     
    if ~strcmp(last_rec_name,origRec_fName) || split_file_flag
        
        
        if strcmp(expType,'operant')
            
            wav_f_num = str2double(origRec_fName_split{6});
            dataFull = audioread(origRec_fName);
            [b,a] = butter(4,500/(fs/2),'high');
            dataFull = filtfilt(b,a,dataFull);
            
        else
            wav_f_num = str2double(origRec_fName_split{2});
            
            if length(origRec_fName_split) > 4
                sample_offset = audio2nlg.total_samples_by_file(wav_f_num) - s.callpos(1);
                data_full_append = s.cut(sample_offset:end);
                dataFull = vertcat(audioread(origRec_fName),data_full_append);
                s.callpos = [s.callpos(1) length(dataFull)];
            else
                dataFull = audioread(origRec_fName);
            end
        end
        data_full_ds = downsample(dataFull,ds_factor); 
    end
    
    callNum = strsplit(callFiles(c).name,{'_','.'});
    callNum = str2double(callNum{end-1});
    
    if callNum == 0 || ~strcmp(last_rec_name,origRec_fName)
        cla(subplot_handles{end})
        plot(subplot_handles{end},(1:length(data_full_ds))/(avi_fs/ds_factor),data_full_ds,'k');
    end
    
    last_rec_name = origRec_fName;
    
    plot(subplot_handles{end},(s.callpos(end,1)+(0:length(data)-1))/avi_fs,data);
    
    if length(data_full_ds)/(avi_fs/ds_factor) > orig_rec_plot_win
       xlim(subplot_handles{end},[s.callpos(end,1)/avi_fs - orig_rec_plot_win/2 s.callpos(end,1)/avi_fs + orig_rec_plot_win/2])
       dataWin = round([s.callpos(end,1) - avi_fs*orig_rec_plot_win/2 s.callpos(end,1) + avi_fs*orig_rec_plot_win/2]);
       dataWin(1) = max(dataWin(1),1);
       dataWin(2) = min(dataWin(2),length(dataFull));
       dataRange = 1.1*[min(dataFull(dataWin(1):dataWin(2))), max(dataFull(dataWin(1):dataWin(2)))];
       ylim(subplot_handles{end},dataRange)
    end
    
    corrected_callpos = (1e3*(s.callpos + cum_samples(wav_f_num))/avi_fs);
    corrected_callpos = avi2nlg_time(audio2nlg,corrected_callpos);
    
    call_ts_nlg = (corrected_callpos + AL_call_offset + audio2nlg.first_nlg_pulse_time)*1e3;
    
    piezo_chunk = cell(1,nLogger);
    logger_k = 1;
    
    call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
        tsData(logger_k).Indices_of_first_and_last_samples(:,1)',...
        tsData(logger_k).Timestamps_of_first_samples_usec,...
        tsData(logger_k).Samples_per_channel_per_file,...
        1e6/nanmean(tsData(logger_k).Estimated_channelFS_Transceiver));
    if isempty(call_ts_nlg_sample)
        break
    end
    
    for logger_k = 1:nLogger
        cla(subplot_handles{logger_k});
        call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
            tsData(logger_k).Indices_of_first_and_last_samples(:,1)',...
            tsData(logger_k).Timestamps_of_first_samples_usec,...
            tsData(logger_k).Samples_per_channel_per_file,...
            1e6/nanmean(tsData(logger_k).Estimated_channelFS_Transceiver));
        if isempty(call_ts_nlg_sample)
            break
        end
        piezo_chunk_idx = call_ts_nlg_sample(1):call_ts_nlg_sample(end);
        try
            piezo_chunk{logger_k} = double(tsData(logger_k).AD_count_int16(piezo_chunk_idx));
            plot_call_spectrogram(piezo_chunk{logger_k},subplot_handles{logger_k},specParams);
        end
        callBounds = 0.5*1e3*[AL_call_offset_length AL_call_offset_length]+[0 abs(diff(corrected_callpos))];
        plot(subplot_handles{logger_k},repmat(callBounds,2,1),specParams.spec_ylims,'k--','LineWidth',4)
        xlim(subplot_handles{logger_k},[0 diff(corrected_callpos + AL_call_offset)]);
    end
    %%
    display(callFiles(c).name);
    
    repeat = 1;
    repeat_k = 1;
    while repeat
        loggerNum = input('call?','s');
        
        switch loggerNum
            
            case 'pause'
                keyboard
            case 'stop'
                fNum = callNums(c);
                save(fullfile(callDir,'current_classify_file_number.mat'),'fNum');
                return
            case 'backup'
                nBackup = input('Backup by how many calls?');
                fNum = callNums(c-nBackup);
                save(fullfile(callDir,'current_classify_file_number.mat'),'fNum');
                return
            otherwise
                loggerNum = str2double(strsplit(loggerNum,','));
                if ~isempty(loggerNum) && all(loggerNum>=1) && all(loggerNum<=nLogger)
                    repeat = 0;
                    noiseStatus = false;
                    if length(loggerNum) == 1
                        manual_al_classify_batNum{c} = tsData(loggerNum).Bat_id;
                    else
                        manual_al_classify_batNum{c} = {tsData(loggerNum).Bat_id};
                    end
                    save_bat_num(s,noiseStatus,manual_al_classify_batNum,manual_bat_ID_fName,call_file_fName)

                elseif loggerNum == 0
                    repeat = 0;
                    noiseStatus = true;
                    manual_al_classify_batNum{c} = 'noise';
                    save_bat_num(s,noiseStatus,manual_al_classify_batNum,manual_bat_ID_fName,call_file_fName)
                    
                elseif loggerNum == -1
                    repeat = 0;
                    noiseStatus = false;
                    manual_al_classify_batNum{c} = 'unidentified';
                    save_bat_num(s,noiseStatus,manual_al_classify_batNum,manual_bat_ID_fName,call_file_fName)
                elseif loggerNum == -2
                    repeat = 0;
                else
                    pause(0.1);
                    if repeat_k < 3
                        sound(data,fs/(2*repeat_k));
                    else
                        startIdx = max(1,s.callpos(end,1) - (orig_rec_plot_win/2)*avi_fs);
                        endIdx = min(length(dataFull),s.callpos(end,1) + (orig_rec_plot_win/2)*avi_fs);
                        sound(dataFull(startIdx:endIdx),fs);
                        repeat_k = 1;
                    end
                    repeat_k = repeat_k + 1;
                end
        end
    end
end


end

function save_bat_num(s,noiseStatus,manual_al_classify_batNum,manual_bat_ID_fName,call_file_fName) %#ok<INUSL>

s.noise = noiseStatus;
save(call_file_fName,'-struct','s');
save(manual_bat_ID_fName,'manual_al_classify_batNum');

end