include("../src/utils.jl")
include("../src/plots.jl")

st = deserialize.([
    "DATA/saved_models/2025-04-10T17:08:28.975new_decay/CuDevice(1)_epoch_008.flux.stats",
    "DATA/saved_models/2025-04-10T21:51:46.666bigger_conv_100/CuDevice(0)_epoch_010.flux.stats",
    "DATA/saved_models/2025-04-11T19:28:37.608skip_chunks100/CuDevice(1)_epoch_010.flux.stats",
    "DATA/saved_models/2025-04-12T17:57:02.547skip_chunks_98/CuDevice(2)_epoch_010.flux.stats",
    "DATA/saved_models/2025-04-12T23:05:36.277pad=125_skip=98/CuDevice(0)_epoch_015.flux.stats",
    "DATA/saved_models/2025-04-13T17:13:17.366pad=125_skip=98_conv=5_dropout=0.15/CuDevice(1)_epoch_015.flux.stats"
]);

plot_mismatch_histograms(st[1].gt_s2s_ranges, st[1].pred_s2s_ranges, st[1].fp_shifts)

st1 = deserialize("DATA/saved_models/2025-04-24T22:47:56.311SOTA_lambda=0.1/CuDevice(0)_epoch_015.flux");
st2 = deserialize("DATA/saved_models/2025-04-24T22:48:21.684SOTA_lambda=0.25/CuDevice(1)_epoch_015.flux");
st3 = deserialize("DATA/saved_models/2025-04-24T22:49:18.431SOTA_lambda=0.5/CuDevice(2)_epoch_015.flux");

p = plot_mismatch_histograms(st3.cds_ranges_true, st3.cds_ranges_model, st3.fp_shifts)