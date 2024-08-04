function plot_temp(df::DataFrame, t::Vector{Float64}; i::Int64, title = "Temperature")
	scatter(
		df.x, df.y, marker_z = t,
		ms = 3, markerstrokewidth = 0,
		markershape = :hexagon,
		label = "", title = "$title, Year: $i",
		color = :viridis, clim = (20, 140))
end


function save_gif(hist, filename = "figs/geothermal.gif", t0 = 2024, fps = 4)
	anim = @animate for (i, step) in enumerate(hist)
		plt = plot_temp(mdp.df, step.s.t, i = t0 + i)
		# savefig(plt, "frame_$(lpad(i, 3, '0')).png")
	end

	return gif(anim, filename, fps = fps)
end
