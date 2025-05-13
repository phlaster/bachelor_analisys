julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.0_pos" -D 2 -b 0.0 -S "pos" &
julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.0_neg" -D 1 -b 0.0 -S "neg" &
julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.0_both" -D 0 -b 0.0 -S "both"

julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.001_pos" -D 2 -b 0.001 -S "pos" &
julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.001_neg" -D 1 -b 0.001 -S "neg" &
julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.001_both" -D 0 -b 0.001 -S "both"

julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.01_pos" -D 2 -b 0.01 -S "pos" &
julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.01_neg" -D 1 -b 0.01 -S "neg" &
julia -t 4 scripts/train_new_model.jl -o DATA/saved_models -s "new_loss_lambda=0.01_both" -D 0 -b 0.01 -S "both"