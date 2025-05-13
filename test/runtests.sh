julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:03.471sota_again_lam=0_pos_strnd/CuDevice(0)_epoch_010.flux" -D 2 -s "pos" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:03.471sota_again_lam=0_pos_strnd/CuDevice(0)_epoch_010.flux" -D 1 -s "neg" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:03.471sota_again_lam=0_pos_strnd/CuDevice(0)_epoch_010.flux" -D 0 -s "both"
sleep 15

julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:03.987sota_again_lam=0_neg_strnd/CuDevice(1)_epoch_010.flux" -D 2 -s "pos" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:03.987sota_again_lam=0_neg_strnd/CuDevice(1)_epoch_010.flux" -D 1 -s "neg" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:03.987sota_again_lam=0_neg_strnd/CuDevice(1)_epoch_010.flux" -D 0 -s "both"
sleep 15

julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:04.035sota_again_lam=0_both_strnd/CuDevice(2)_epoch_010.flux" -D 2 -s "pos" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:04.035sota_again_lam=0_both_strnd/CuDevice(2)_epoch_010.flux" -D 1 -s "neg" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T00:09:04.035sota_again_lam=0_both_strnd/CuDevice(2)_epoch_010.flux" -D 0 -s "both"
sleep 15

julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:56.426sota_again_lam=0.2_pos_strnd/CuDevice(0)_epoch_010.flux" -D 2 -s "pos" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:56.426sota_again_lam=0.2_pos_strnd/CuDevice(0)_epoch_010.flux" -D 1 -s "neg" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:56.426sota_again_lam=0.2_pos_strnd/CuDevice(0)_epoch_010.flux" -D 0 -s "both"
sleep 15

julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:56.131sota_again_lam=0.2_neg_strnd/CuDevice(1)_epoch_010.flux" -D 2 -s "pos" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:56.131sota_again_lam=0.2_neg_strnd/CuDevice(1)_epoch_010.flux" -D 1 -s "neg" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:56.131sota_again_lam=0.2_neg_strnd/CuDevice(1)_epoch_010.flux" -D 0 -s "both"

julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:55.828sota_again_lam=0.2_both_strnd/CuDevice(2)_epoch_010.flux" -D 2 -s "pos" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:55.828sota_again_lam=0.2_both_strnd/CuDevice(2)_epoch_010.flux" -D 1 -s "neg" &
julia -t 4 scripts/evaluate_model.jl -d "DATA/saved_models/2025-05-06T07:25:55.828sota_again_lam=0.2_both_strnd/CuDevice(2)_epoch_010.flux" -D 0 -s "both"