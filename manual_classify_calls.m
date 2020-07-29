function manual_classify_calls(wd,call_str)

switch call_str
    case 'Call'
        callDir = fullfile(wd,'Analyzed_auto');
    case 'Echo'
        callDir = fullfile(wd,'Analyzed_auto_echo');
end


specWindow = kaiser(512,0.5);
specOverlap = 450;
freqPoints = 2^10;
orig_rec_plot_win = 5;
recVar = 'cut';
recVar_full = 'recsGroup';

wav_mat_file = input('wav (1) or mat (2) file?');

callFiles = dir(fullfile(callDir,['*' call_str '*.mat']));
callNums = 1:length(callFiles);
call_file_name_length = 21;

if exist(fullfile(callDir,'current_classify_file_number.mat'),'file')
    f = load(fullfile(callDir,'current_classify_file_number.mat'));
    fNum = find(callNums == abs(min(callNums-f.fNum)));
else
    fNum = 1;
end

replotFlag = true;
    
nCalls = length(callFiles);
for c = fNum:nCalls
    s = load(fullfile(callDir,callFiles(c).name));
    data = s.(recVar);
    fs = min(s.fs,200e3);
    sound(data,fs);
    origRec_fName = strsplit(callFiles(c).name,'_');
    origRec_file_name = strjoin(origRec_fName(1:end-2),'_');
    if length(origRec_file_name) > call_file_name_length
        origRec_file_name = origRec_file_name(1:call_file_name_length);
    end
    if wav_mat_file == 1
        origRec_fName = fullfile(wd, [origRec_file_name '.WAV']);
        dataFull = audioread(origRec_fName);
    elseif wav_mat_file == 2
        origRec_fName = fullfile(wd, [origRec_file_name '.mat']);
        d = load(origRec_fName);
        dataFull = d.(recVar_full);
    end
    
    callNum = strsplit(callFiles(c).name,{'_','.'});
    callNum = str2double(callNum{end-1});
    
    subplot(2,1,1)
    cla
    spectrogram(data,specWindow,specOverlap,freqPoints,s.fs,'yaxis');
    subplot(2,1,2);
    hold on
    if callNum == 0 || replotFlag
        cla
        plot((1:length(dataFull))/s.fs,dataFull,'k');
        replotFlag = false;
    end
    plot((s.callpos(end,1)+(0:length(data)-1))/s.fs,data);
    
    if length(dataFull)/s.fs > orig_rec_plot_win
       xlim([s.callpos(end,1)/s.fs - orig_rec_plot_win/2 s.callpos(end,1)/s.fs + orig_rec_plot_win/2])
    end
    
    display(callFiles(c).name);
  
    repeat = 1;
    repeat_k = 1;
    while repeat
        class = input('call?','s');
        switch class
            case '1'
                repeat = 0;
                s.noise = false;
                save(fullfile(callDir,callFiles(c).name),'-struct','s');
            case '0'
                repeat = 0;
                s.noise = true;
                save(fullfile(callDir,callFiles(c).name),'-struct','s');
            case 'stop'
                fNum = callNums(c);
                save(fullfile(callDir,'current_classify_file_number.mat'),'fNum');
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

end