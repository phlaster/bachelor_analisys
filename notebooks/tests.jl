### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 4647ed1e-e913-4874-b650-603d7c6d868c
begin
	using Pkg
	Pkg.activate("../")
end

# ╔═╡ bd5e0cbf-f490-42e9-8f59-dadfd0ae6ec8
using SimplePlutoInclude, ProgressLogging

# ╔═╡ db786127-3961-4429-b7d2-465132719c5d
# ╠═╡ show_logs = false
using Flux, Statistics, ProgressMeter, LuxCUDA

# ╔═╡ 3bc7cb9c-6b0b-4101-9a29-a814ffa768bf
@plutoinclude "../src/utils.jl"

# ╔═╡ f0760cd9-7f5a-4f5e-8878-83e891b51dc1
dirs = "../DATA/genomes/genomes/" .* readdir("../DATA/genomes/genomes")[1]

# ╔═╡ c026b3f6-6e7f-425c-b984-bfd13e4f86d9
pad = 50

# ╔═╡ 6e1f0a72-921f-40e6-807e-6c5e939f587e
window_size = 2 * pad + 1

# ╔═╡ 733b288f-aec7-4af8-95f2-65d921d5a0d5
gd = GenomeDataset([dirs]; pad=pad)

# ╔═╡ cfb25cf7-fe96-4d7e-9fea-6bd6ae1c4679
prepared_padded_seq, prepared_labels = gd[1]

# ╔═╡ c93d0f8b-0af8-4fc5-a103-25c9de369a3d
begin
	x = permutedims(prepared_padded_seq, (2, 1))  # 12×5
	x = reshape(x, size(x, 1), size(x, 2), 1)    # 12×5×1
	y = prepared_labels  # 4×6 (classes × positions)
end

# ╔═╡ 59f30bcd-e948-483a-bd79-25168e7ad21d
device = gpu_device()

# ╔═╡ 8cdcc783-fac9-4346-af18-44a6767a7853
model = Chain(
    Conv((window_size,), 5 => 16, pad=0, relu),
    Dropout(0.2),
    Conv((3,), 16 => 32, pad=1, relu),
    Dropout(0.2),
    Conv((1,), 32 => 4)
) |> device

# ╔═╡ 996c558b-1df5-4d70-b946-83dac561e0d2
begin
	opt_state = Flux.setup(Adam(0.01), model)
	losses = []
	@progress for epoch in 1:10
	    loss, grads = Flux.withgradient(model) do m
	        y_hat = m(x |> device)
	        Flux.logitcrossentropy(permutedims(dropdims(y_hat; dims=3), (2, 1)), y |> device)
	    end
	    Flux.update!(opt_state, model, grads[1])
	    push!(losses, loss)
	end
end

# ╔═╡ c6e0b4aa-0bd7-4dae-9cf9-a774108a7895
y_pred = model(x |> device)

# ╔═╡ b811e793-9c8a-4a46-a351-04a7fcee0fbc
probs = softmax(dropdims(y_pred; dims=3), dims=2)

# ╔═╡ adee8565-ce87-4b9f-9eac-3b87d87e4719
@time pred_classes = argmax.(eachrow(probs |> cpu)) .- 1

# ╔═╡ 2191edd3-fb23-4b47-8cb7-feaa92590aa0
accuracy = mean(pred_classes .== Flux.onecold(prepared_labels, 0:3))

# ╔═╡ Cell order:
# ╠═bd5e0cbf-f490-42e9-8f59-dadfd0ae6ec8
# ╠═4647ed1e-e913-4874-b650-603d7c6d868c
# ╠═3bc7cb9c-6b0b-4101-9a29-a814ffa768bf
# ╠═db786127-3961-4429-b7d2-465132719c5d
# ╠═f0760cd9-7f5a-4f5e-8878-83e891b51dc1
# ╠═c026b3f6-6e7f-425c-b984-bfd13e4f86d9
# ╠═6e1f0a72-921f-40e6-807e-6c5e939f587e
# ╠═733b288f-aec7-4af8-95f2-65d921d5a0d5
# ╠═cfb25cf7-fe96-4d7e-9fea-6bd6ae1c4679
# ╠═c93d0f8b-0af8-4fc5-a103-25c9de369a3d
# ╠═59f30bcd-e948-483a-bd79-25168e7ad21d
# ╠═8cdcc783-fac9-4346-af18-44a6767a7853
# ╠═996c558b-1df5-4d70-b946-83dac561e0d2
# ╠═c6e0b4aa-0bd7-4dae-9cf9-a774108a7895
# ╠═b811e793-9c8a-4a46-a351-04a7fcee0fbc
# ╠═adee8565-ce87-4b9f-9eac-3b87d87e4719
# ╠═2191edd3-fb23-4b47-8cb7-feaa92590aa0
