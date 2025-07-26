using KomaMRI, MAT

## CREATE PHANTOM

function mm_to_idx(x_mm)
    return round(Int, (x_mm + phantom_radius_mm) / voxel_size) + 1
end

function idx_to_mm(i::Int)
    return (i - 1) * voxel_size - phantom_radius_mm
end

voxel_size = 5e-3
phantom_radius_mm = 60.0e-3       # physical phantom = 200 mm diameter
grid_mm = -phantom_radius_mm:voxel_size:phantom_radius_mm 
phantom_size = (length(grid_mm), length(grid_mm), length(grid_mm))

num_spheres = 14
vial_radius_mm = 7.5e-3  # 15 mm diameter

T1_vals = [2640, 2292, 1923, 1489, 1245, 1004, 733.9,
           533.1, 400.3, 261.0, 189.8, 154.7, 102.1, 79.65]
T2_vals = [1044, 623.9, 428.3, 258.4, 186.1, 137.0, 89.52,
           62.82, 43.84, 27.28, 19.24, 15.44, 10.05, 7.79]

T1_vals = T1_vals/1000   # convert to s
T2_vals = T2_vals/1000

# vial centers
outer_radius_mm = 50.0e-3
outer_centers = [
    (
        0.0,
        outer_radius_mm * cos(2π * (i-1) / 10),
        outer_radius_mm * sin(2π * (i-1) / 10)
    ) for i in 1:10
]

inner_offset = 20.0e-3
inner_centers = [
    ( 0.0, -inner_offset, -inner_offset),
    ( 0.0, inner_offset, -inner_offset),
    ( 0.0, inner_offset,  inner_offset),
    (0.0, -inner_offset,  inner_offset)
]

centers_mm = vcat(outer_centers, inner_centers)

T1_vol = zeros(Float64, phantom_size)
T2_vol = zeros(Float64, phantom_size)

for i in 1:num_spheres
    cx, cy, cz = centers_mm[i]

    xi_min = clamp(mm_to_idx(cx - vial_radius_mm), 1, phantom_size[1])
    xi_max = clamp(mm_to_idx(cx + vial_radius_mm), 1, phantom_size[1])
    yi_min = clamp(mm_to_idx(cy - vial_radius_mm), 1, phantom_size[2])
    yi_max = clamp(mm_to_idx(cy + vial_radius_mm), 1, phantom_size[2])
    zi_min = clamp(mm_to_idx(cz - vial_radius_mm), 1, phantom_size[3])
    zi_max = clamp(mm_to_idx(cz + vial_radius_mm), 1, phantom_size[3])

    for xi in xi_min:xi_max, yi in yi_min:yi_max, zi in zi_min:zi_max
        x_mm = idx_to_mm(xi)
        y_mm = idx_to_mm(yi)
        z_mm = idx_to_mm(zi)
        r2 = (x_mm - cx)^2 + (y_mm - cy)^2 + (z_mm - cz)^2
        if r2 <= vial_radius_mm^2
            T1_vol[xi, yi, zi] = T1_vals[i]
            T2_vol[xi, yi, zi] = T2_vals[i]
        end
    end
end

inds = findall(T1_vol .> 0)

x = Float64[]
y = Float64[]
z = Float64[]
T1 = Float64[]
T2 = Float64[]

for idx in inds
    i, j, k = Tuple(idx)
    push!(x, grid_mm[i])
    push!(y, grid_mm[j])
    push!(z, grid_mm[k])
    push!(T1, T1_vol[i, j, k])
    push!(T2, T2_vol[i, j, k])   
end

phantom = Phantom{Float64}(
    x = x,
    y = y,
    z = z,
    T1 = T1,
    T2 = T2,
)

p = plot_phantom_map(phantom, :T1)
display(p)

# save phantom
phantom_dict = Dict(
    "x" => collect(x),  
    "y" => collect(y),
    "z" => collect(z),
    "T1" => collect(T1),
    "T2" => collect(T2),
)

MAT.matwrite("vialphantom3D_0.5mm_voxel.mat", Dict("phantom" => phantom_dict))


## SIMULATE 

sys = Scanner()
seq = read_seq("sequences/mpf_001_PhantomStudy_short.seq")

sim_params = KomaMRICore.default_sim_params()
sim_params["return_type"] = "raw"

raw = simulate(phantom, seq, sys; sim_params=sim_params)

p = plot_signal(raw)
display(p)

signal = Array{ComplexF64}(undef, 248, 1, 1000)

for i in 1:1000
    signal[:, 1, i] = raw.profiles[i].data
end

matwrite("komamri_signal.mat", Dict("signal" => signal))


## TRAJECTORY

_, ktraj = get_kspace(seq)

ktraj_reshaped = reshape(ktraj, 248, 1000, 3)  # (248, 1000, 3)
traj = permutedims(ktraj_reshaped, (3, 1, 2))  # (3, 248, 1000)

matwrite("komamri_traj.mat", Dict("traj" => traj))

p = plot_kspace(seq)
display(p)
