using Plots
anim = @animate for i in 1:50
    plot(sin, 0, i * 2pi / 10)
end when i > 30
gif(anim, "movie.gif",fps=20)