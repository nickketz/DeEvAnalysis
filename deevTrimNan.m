function deevTrimNan(dtype,subs,adFile)
% step through each subject in exper, and trim nans from dtype data
% kind of a hack to fix the fact that process_ft_data doesn't take
% time ranges per condition.  This function goes through and
% trims the data of surrounding nans in the temporal dimension
%
% input:
%   dtype = str; data type to process (i.e. pow, fourier, wpli etc)
%   subs = cell; sub strings to process, if empty do all in exper
%   adFile = str; analysis details file def:'/home/nike3851/work/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'
%
%

if ~exist('adFile','var')   adFile = '/home/nike3851/work/DeEv_EEG/analysis/data/DEEV/EEG/Sessions/long/ft_data/rFix_OL_CL_rlOL_rlCL_Enc_eq0_art_ftAuto/tfr_wavelet_w4_fourier_3_50/analysisDetails.mat'; end

load(adFile);

if ~exist('subs','var')     subs = exper.subjects;          end
if ~iscell(subs)            subs = {subs};                  end

switch dtype
    case {'fourier'}
        param = 'fourierspctrm';
    case {'pow'}
        param = 'powspctrm';
    otherwise
        error('\nunrecognized data type: %s',dtype);
end



for isub = 1:length(subs)
    fprintf('\nSubject %s:',subs{isub});
    for ises = 1:length(exper.sesStr)
        try
            pth = fullfile(dirs.saveDirProc,subs{isub},exper.sesStr{ises});
            myfiles = dir(fullfile(pth,['*' dtype '*.mat']));
            if isempty(myfiles)
                warning(sprintf('\nno %s files found for %s, %s',dtype,subs{isub},exper.sesStr{ises}));
                continue;
            end
            myfiles = {myfiles.name};
            for ifile = 1:length(myfiles)
                fprintf('\n  loading %s...',myfiles{ifile});
                load(fullfile(pth,myfiles{ifile}));
                ndim = length(size(freq.(param)));
                %assume time in last dim
                nanind = isnan(freq.(param));
                for idim = 1:ndim-1, nanind = squeeze(sum(nanind,1)); end
                nanind = nanind>0;
                
                if sum(~nanind)==0 fprintf('no non-nan time points, skipping'); continue; end
                if sum(~nanind)==length(nanind) fprintf('all time points are non-nan, skipping'); continue; end
                
                tgood = freq.time(~nanind);
                fprintf('trimming %.02fsto%.02fs...',min(tgood),max(tgood));
                evalc('freq = ft_selectdata(struct(''latency'',[min(tgood) max(tgood)]),freq);');
                %save trimmed data
                fprintf('saving...');
                save(fullfile(pth,myfiles{ifile}),'freq');
                fprintf('done');
            end
        catch ME
            fprintf('%s\n',ME.message);
            fprintf('error! skipping file');
        end
    end
end