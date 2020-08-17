function plot_call_spectrogram(callWF,axisHandle,specParams)
if ~isfield(specParams,'spec_caxis_factor')
    specParams.spec_caxis_factor = 0.75;
end

if length(callWF) <= specParams.spec_win_size
    return
end

[~,f,t,ps] = spectrogram(callWF,gausswin(specParams.spec_win_size),specParams.spec_overlap_size,specParams.spec_nfft,specParams.fs,'yaxis');
if length(t) <= 1
    return
end
f = f*1e-3; % frequency in kHz
t = t*1e3;
t = [0 t(end-1)];
ps = 10*log10(ps);
ps(isinf(ps)) = NaN;
imagesc(axisHandle,t,f,ps)
cmap = spec_cmap();
colormap(cmap);
ylim(axisHandle,specParams.spec_ylims);
caxis(axisHandle,[min(ps(:))*specParams.spec_caxis_factor max(ps(:))]);
set(axisHandle,'YDir','normal')

end