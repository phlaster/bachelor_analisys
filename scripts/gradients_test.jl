using Flux
using Flux.Losses: binary_focal_loss
using Zygote
using LuxCUDA

model = Chain(Dense(10, 10), sigmoid) |> gpu

x = randn(Float32, 10) |> gpu
y = rand(Bool, 10) |> gpu


function spatial_penalty(y; window=4)
    N = length(y)
    d_max = min(window, N - 1)
    sum(@inbounds sum(@view(y[1:N-d]) .* @view(y[d+1:N])) for d in 1:d_max)
end

function floss_with_spatial_penalty(ŷ, y; lambda=1.f0, window=4)
    fl = binary_focal_loss(ŷ, y)
    sp = spatial_penalty(ŷ, window=window)
    return fl + lambda * sp
end

loss0, gs_without_penalty = withgradient(model) do m
    ŷ = m(x)
    binary_focal_loss(ŷ, y)
end

loss1, gs_zero_lambda = withgradient(model) do m
    ŷ = m(x)
    floss_with_spatial_penalty(ŷ, y, lambda=0.f0)
end

loss2, gs_with_penalty = withgradient(model) do m
    ŷ = m(x)
    floss_with_spatial_penalty(ŷ, y)
end

(loss0, gs_without_penalty[1].layers[1].weight[1:5]) |> println
(loss1, gs_zero_lambda[1].layers[1].weight[1:5]) |> println
(loss2, gs_with_penalty[1].layers[1].weight[1:5]) |> println

# Check if gradients differ
are_equal = gs_with_penalty.grad[1].layers[1].weight ≈ gs_without_penalty.grad[1].layers[1].weight
println("Gradients are equal? $are_equal")

