function [logger_delay,callPos,logger_xcorr] = get_audio_AL_delay(cut_call_data,tsData,audio2nlg,manual_al_classify_batNum)

al_fs = round(1/(1e-6*unique([tsData.Sampling_period_usec_Logger])));
wav_fs = unique([cut_call_data.fs]);

for call_k = 1:length(cut_call_data)
    cut_call_data(call_k).batNum = manual_al_classify_batNum{call_k};
    cut_call_data(call_k).noise = strcmp(manual_al_classify_batNum{call_k},'noise');
end

max_delay_idx = ~strcmp(manual_al_classify_batNum,'noise');
cut_call_data = cut_call_data(max_delay_idx);
nCalls = length(cut_call_data);

batNums = {tsData.Bat_id};

assert(length(al_fs) == 1 && length(wav_fs) == 1)

[resample_N,resample_D] = rat(wav_fs/al_fs);
envelope_window_s = 5e-3;

AL_call_offset_length = 0.1;
AL_call_offset = 1e3*0.5*AL_call_offset_length*[-1 1];
f_bounds = [1 3]*1e3;
[b,a]=butter(3,f_bounds/(al_fs/2),'bandpass');

callPos = vertcat(cut_call_data.corrected_callpos);
logger_delay = nan(1,nCalls);
logger_xcorr = nan(1,nCalls);

for c = 1:nCalls
    
    call_wav_env = padarray(cut_call_data(c).cut,0.5*AL_call_offset_length*wav_fs,0,'both');
    call_wav_env = zscore(resample(envelope(call_wav_env,envelope_window_s*wav_fs,'rms'),resample_D,resample_N));
    call_wav_env_L = length(call_wav_env);
    
    corrected_callpos = callPos(c,:);
    call_ts_nlg = (corrected_callpos + AL_call_offset + audio2nlg.first_nlg_pulse_time)*1e3;
    
    if ischar(cut_call_data(c).batNum)
        logger_k = find(strcmp(batNums,cut_call_data(c).batNum));
    else
        continue
    end
    
    call_ts_nlg_sample = get_voltage_samples_for_Nlg_timestamps(call_ts_nlg(:),...
        tsData(logger_k).Indices_of_first_and_last_samples(:,1)',...
        tsData(logger_k).Timestamps_of_first_samples_usec,...
        tsData(logger_k).Samples_per_channel_per_file,...
        1e6/nanmean(tsData(logger_k).Estimated_channelFS_Transceiver));
    if isempty(call_ts_nlg_sample) || any(isnan(call_ts_nlg_sample))
        continue
    end
    piezo_chunk_idx = call_ts_nlg_sample(1):call_ts_nlg_sample(end);
    call_logger_chunk = double(tsData(logger_k).AD_count_int16(piezo_chunk_idx));
    
    call_al_env = zscore(envelope(filtfilt(b,a,call_logger_chunk),envelope_window_s*al_fs,'rms'));
    call_al_env_L = length(call_al_env);
    N = min(call_al_env_L,call_wav_env_L);
    
    r = xcorr(call_al_env(1:N),call_wav_env(1:N));
    [logger_xcorr(c),max_delay_idx] = max(r);
    logger_delay(c) = (max_delay_idx - N)/al_fs;
    
%     plot(call_al_env(1:N))
%     plot(call_wav_env(1:N))
%     input('?')
%     cla
end

