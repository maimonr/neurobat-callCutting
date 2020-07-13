function findcalls_session(wd,fs,varargin)
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

pnames = {'fileType', 'outputDir', 'dataVarName','findcall_inputs','audio_file_filter','filter_raw_data','expType','clippingCorrection'};
dflts  = {'wav', fullfile(wd, 'Analyzed_auto'),'recsGroup',{},'*',false,'communication',false};
[fileType,outputDir,dataVarName,findcall_inputs,audio_file_filter,filter_raw_data,expType,clippingCorrection] = internal.stats.parseArgs(pnames,dflts,varargin{:});

switch expType
    case 'communication'
        rec_files = dir(fullfile(wd,[audio_file_filter '.' fileType]));
    case 'operant'
        rec_files = dir(fullfile(wd,[audio_file_filter '.' fileType]));
        rec_file_numbers = cellfun(@(x) str2double(x(strfind(x,'mic1_')+length('mic1_'):strfind(x,'.wav')-1)),{rec_files.name});
        [~,wav_file_order] = sort(rec_file_numbers,'ascend');
        rec_files = rec_files(wav_file_order);
end

debugFlag = false;
if any(strcmp(findcall_inputs,'debug'))
    if findcall_inputs{find(strcmp(findcall_inputs,'debug'))+1}
        debugFlag = true;
    end
end

if ~debugFlag && ~isfolder(outputDir)
    mkdir(outputDir);
end
n_files = length(rec_files);

for fln = 1:n_files
    filename=rec_files(fln).name;
    disp(filename)
    disp(['Analyzing file: ' num2str(fln) ' of ' num2str(n_files)])
    disp('...')
    % load data and calculate envelope
    data_raw = load_audio_data(wd,filename,fileType,dataVarName,filter_raw_data,clippingCorrection);
    filename = filename(1:end-4);
    if fln < n_files
        next_filename = rec_files(fln+1).name;
    else
        next_filename = [];
    end
    findcall_inputs = {findcall_inputs{:} 'filename' filename 'next_filename' next_filename};
    findcalls(data_raw,wd,fs,findcall_inputs{:})
end

end



function data_raw = load_audio_data(wd,filename,fileType,dataVarName,filterFlag,clippingCorrection)

switch fileType
    case 'wav'
        [data_raw,fs] = audioread(fullfile(wd, filename));
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

if clippingCorrection
   clippingIdx = find(diff(data_raw) < -0.9) + 1;
   while ~isempty(clippingIdx)
       data_raw(clippingIdx) = data_raw(clippingIdx) + 1;
       clippingIdx = find(diff(data_raw) < -0.9) + 1;
   end
end

if filterFlag
    [b,a] = butter(4,500/(fs/2),'high');
    data_raw = filtfilt(b,a,data_raw);
end

end