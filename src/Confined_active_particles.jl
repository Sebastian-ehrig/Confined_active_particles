using NeuralPDE, Lux, Optimization, OptimizationOptimJL
import ModelingToolkit: Interval
using ModelingToolkit
using IntervalSets
using Makie  # great package
using GLMakie  # great package
using MeshIO  # good package
using FileIO  # good package
using Meshes
using GeometryBasics  # great package
using Statistics
using LinearAlgebra


GLMakie.activate!()
GLMakie.set_window_config!(
    framerate = 10,
    title = "Confined active particles"
)


UpFolder = pwd();

# import an STL mesh, returning a PATCH-compatible face-vertex structure
namestructure = "ellipsoid_x4"
mesh_loaded = FileIO.load("assets/ellipsoid_x4.stl")

##############################################################


# TODO: implement this: https://docs.makie.org/v0.19.1/documentation/nodes/index.html#the_observable_structure
# TODO
"""
    simulate_on_mesh(mesh_surf, form_particle, n_particle)

The goal is to create a package where I can overgive a mesh surface and the particle form and then simply simulate the particle behaviour on the mesh surface.
"""
function simulate_on_mesh(mesh_surf, form_particle, n_particle)

end


# TODO: implement this correctly
# # # % Plot the particles
# using GLMakie, Random
# GLMakie.activate!()
# num_part = 200
# initms = 8 * rand(num_part)
# particles = Makie.Observable(initms) # this is the variable that will change
# # first frame, initial plot
# fig, ax = scatter(2 * rand(num_part), rand(num_part), markersize = particles,
#     color = initms, colormap = (:Greek, 0.75), strokewidth = 0.5,
#     strokecolor = :white,
#     figure = (resolution = (1200, 800), fontsize = 22),
#     axis = (xlabel = "x", ylabel = "y",))
# limits!(ax, 0, 2, 0, 1)
# # the animation is done by updating the Observable values
# record(fig, joinpath("assets", "confined_active_particles.mp4"),
#     framerate = 24, profile = "main") do io
#     for i in 1:0.1:8
#         particles[] = i * initms
#         recordframe!(io)  # record a new frame
#     end
# end


"""
    clip11(x)

"""
function clip11(x)
    if x < -1
        return -1
    elseif x > 1
        return 1
    else
        return x
    end
end


function repel(particles_node, N)
    particles = particles_node[]
    @inbounds for i in 1:N
        ftot = Vec3f0(0)
        p1 = particles[i]
        for j in 1:N
            if i != j
                p2 = particles[j]
                Δσ = acos(clip11(dot(p1, p2))) # great circle distance
                ftot += (p1 - p2)/max(1e-3, Δσ^2)
            end
        end
        particles[i] = normalize(p1 + 0.001 * ftot)
    end
    particles_node[] = particles
end


function addparticle!(particles, colors, nparticles)
    # ! You change the value of an Observable with empy index notation:
    nparticles[] = nparticles[] + 1
    particles[][nparticles[]] = normalize(randn(Point3f0))
    colors[][nparticles[]] = to_color(:green)
    particles[] = particles[]
    colors[] = colors[]
end

scene = Makie.Scene(resolution = (400,400), show_axis = false);
mesh!(scene, Sphere(Point3f0(0), 1f0), color = :gray)

#f, ax, pl = Makie.mesh(mesh_loaded, axis=(type=Axis3,))  # plot the mesh

max_particles = 5000
# An `Observable` is a mutable container of an object of type `T`.
# `T` can be any type.
particles = Makie.Observable(fill(Point3f0(NaN), max_particles))
colors = Makie.Observable(fill(RGBAf(0, 0, 0, 0), max_particles))
meshscatter!(scene, particles, color = colors, markersize = 0.05)
nparticles = Makie.Observable(0)
for i=1:10
    addparticle!(particles, colors, nparticles)
end
update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))
# scene.center = false # don't reset the camera by display
N = 1000
record(scene, "assets/output.mp4", 1:N) do iter
    isodd(iter) && addparticle!(particles, colors, nparticles)
    repel(particles, nparticles[])
end







"""
    calculate_order_parameter(v_order, tt, v_tp)

"""
function calculate_order_parameter(v_order, tt, v_tp)
    for i=1:n
        v_norm[i,:]=v_tp[i,:]/norm(v_tp[i,:])
    end
    v_order[tt]=(1/n)*norm(sum(v_norm))

    return v_order
end


"""
    find_nonzero_index(c::Array)

"""
function find_nonzero_index(c::Array)
    a = similar(c, Int)
    count = 1
    @inbounds for i in eachindex(c)
        a[count] = i
        count += (c[i] != zero(eltype(c)))
    end
    return resize!(a, count-1)
end


"""
    dim_data(V::Array{Float32, 2}, F, dim::Int)

"""
function dim_data(V::Array{Float32, 2}, F, dim::Int)
    return reduce(vcat,transpose.([V[F[:,1],dim],V[F[:,2],dim],V[F[:,3],dim]]))'
end


"""
    vec_of_vec_to_array(V::Array{Float32, 2})

transform vector of vectors to matrix
"""
function vec_of_vec_to_array(V)
    reduce(vcat,transpose.(V))
end


# NOTE: maybe this is not correct
"""
    get_index_binary(Dist2faces::Array)

Faces with closest point to r(i,:) and associated normal vectors
"""
function get_index_binary(_Dist2faces::Array)
    return first.(Tuple.(sum(findall(x->x==minimum(_Dist2faces), _Dist2faces), dims=2)))
end


"""
    get_particle(_face_coord_temp::Array, i::Int)

"""
function get_particle(_face_coord_temp::Array, i::Int)
    particle = cat(dims=2,ones(size(_face_coord_temp[:,1,3]))*r[i,1],
        ones(size(_face_coord_temp[:,1,3]))*r[i,2],
        ones(size(_face_coord_temp[:,1,3]))*r[i,3]
        )
    return particle
end


"""
    get_particle_position(_face_coord_temp::Array, N_temp, i::Int)

"""
function get_particle_position(_face_coord_temp::Array, N_temp, i::Int)
    particle = get_particle(_face_coord_temp, i)
    len_face_coord = length(_face_coord_temp[:,1,1])
    reshape_it = reshape(_face_coord_temp[:,1,:],len_face_coord,3)
    placeholder = sum((particle-reshape_it).*N_temp,dims=2)
    p0s = particle - cat(placeholder, placeholder, placeholder, dims=2).*N_temp
    return p0s
end


##############################################################
# Define perpendicular projection functions, adapted from "PhysRevE 91 022306"
##############################################################

# P_perp does a normal projection of the vector b on the plane normal to a
P_perp(a, b) = (b-(sum(b.*a,dims=2)./(sqrt.(sum(a.^2,dims=2)).^2)*ones(1,3)).*a)

# P_plan does a projection of the vector b on vector normal to a1 in the plane normal to a
P_plan(a,b,a1) = ((sum(b.*a,dims=2)./(sqrt.(sum(a.^2,dims=2)))*ones(1,3)).*[
    a[:,2].*a1[:,3]-a[:,3].*a1[:,2],-a[:,1].*a1[:,3]+a[:,3].*a1[:,1],a[:,1].*a1[:,2]-a[:,2].*a1[:,1]
    ])


##############################################################
# PLOTTING DONE BY MAKIE.jl
##############################################################

scene = Makie.Scene(resolution = (400,400));
f, ax, pl = Makie.mesh(mesh_loaded, axis=(type=Axis3,))  # plot the mesh
wireframe!(ax, mesh_loaded, color=(:black, 0.2), linewidth=2, transparency=true)  # only for the asthetic


##############################################################
# Get the data from the mesh
##############################################################

N = mesh_loaded.normals
# faces = GeometryBasics.decompose(TriangleFace{Int}, mesh_loaded)  # return the faces of the mesh
vertices = GeometryBasics.coordinates(mesh_loaded)  # return the vertices of the mesh
faces = reshape(1:length(vertices), 3, :)'  # return the faces of the mesh
vertices = vec_of_vec_to_array(vertices)

Vertex = vertices
Faces = faces


##############################################################
# When we plot surfaces, we are interested only in some part of the
# geometry, so we select the faces which have at least one vertex in the
# area of interest
V_plot = ones(length(Faces[:,1]),1); # Will store the num of the faces in the ROI
face_num_plot=[];


Faces_coord = cat(dim_data(Vertex, Faces, 1), dim_data(Vertex, Faces, 2), dim_data(Vertex, Faces, 3), dims=3)


##############################################################
# Define Model Parameters
##############################################################

const visible_off = 0; # 0 to display figure while plotting them, 1 to just save them without being displayed first

const movie_making = 1; #1 for movies, 0 for saving images

const num_part = 200; # number particles
const v0 = 0.1; # velocity of particles
const v0_next = 0.1; # if velocity is to be changed after a certain number of timesteps

const num_step = 300; # Number of timesteps

const Area=452.389; # total surface area

const sigma_n = 5/12; # particle radius

const phi_pack=(num_part*pi*(sigma_n^2))/Area; # packing-fraction (~density)

# Define parameters for the force calculations:
# ---------------------------------------------
const r_adh = 1; # Cut off radius for adhesive forces
const k = 10.; # Elastic constant for repulsive forces
const k_next = 10.; # value of elastic constant after a certain number of timesteps
const k_adh = 0.75; # Adhesive forces

const mu = 1.; # arbitrary mass of particle
const tau = 1; #time of relaxation of collision between 2 particles

const dt = 0.01*tau; # Step size
const plotstep = 0.1/dt; # Number of calculation between plot

# const a = 0.5; # arbitrary value
const b = v0_next/(2*sigma_n*mu*k_next); # for calculating coupling constant
const J = b/tau; # Coupling constant of orientation vectors
const noise=0.002/tau;  #value of the noise amplitude
const hi=noise/(2*sqrt(dt)); # Amplitude of the noise on orientation

# Parameters for making movies:
#------------------------------
const rho = 6; # arbitrary value for scaling the vector on plotting. To automatize, maybe from file name value
const scale=0.1*rho; # Scale the vector for making movies
const fps = 25; # Frame per second of movie

time_points = Int(num_step/plotstep)
v_order = zeros(time_points, 1)
en_v = zeros(time_points, 1)
tot_F = zeros(time_points, 1)


# particle_info = struct;  # TODO: habe ich rausgenommen, weil ich nicht weiß, was das macht

# name_data = [namestructure,'_N_',num2str(num_part),'_rho_',num2str(rho)...
#         ,'_F_rep_',num2str(k),'_F_adh_',num2str(k_adh),'_tau_',num2str(tau),'_b_',num2str(b),...
#         '_s',num2str(v0),'_t',num2str(num_step),'_fps',num2str(25)];
    
# % Define folder structure and pre-process meshes:
# % -----------------------------------------------
folder_plots = joinpath(UpFolder, "images", namestructure)
folder_particle_simula = joinpath(UpFolder, "simulation_structure")


##############################################################
# Create index matrix of neighbourg faces to each face
##############################################################

F_neighbourg = fill(NaN, length(Faces[:,1]), 100)  # matrix of neighbourg faces

maximum_neighbour = 0  # maximumimum number of neighbourg faces
# % Loop for to search all neighbours for each faces within a radius 
# % "radius_search" centered around the isobarycenter of the face. 
# % The face are considered within the radius if at least one of the
# % vertex is within. radius_search = displacement of particle + distance
# % between isobarcenter and verteces of face considered
for i = 1:length(Faces[:,1])
    center_faces = [Faces_coord[i,1,1]+Faces_coord[i,2,1]+Faces_coord[i,3,1],
        Faces_coord[i,1,2]+Faces_coord[i,2,2]+Faces_coord[i,3,2],
        Faces_coord[i,1,3]+Faces_coord[i,2,3]+Faces_coord[i,3,3]
        ]/3

    extra_dist = sqrt((center_faces[1]-Faces_coord[i,1,1])^2+(center_faces[2]-Faces_coord[i,1,2])^2+(center_faces[3]-Faces_coord[i,1,3])^2)

    radius_search = extra_dist # +dist_motion  # TODO: warum ist in dieser Gleichung dist_motion nicht definiert?
    Faces2center = Faces_coord - cat(center_faces[1]*ones(size(Faces)),
        center_faces[2]*ones(size(Faces)),center_faces[3]*
        ones(size(Faces)), dims=3)

    # % Norm^2 vector all faces verteces to vertex 1 of this face
    Faces2center = Faces2center[:,:,1].*Faces2center[:,:,1] + Faces2center[:,:,2].*Faces2center[:,:,2] + Faces2center[:,:,3].*Faces2center[:,:,3]
    # assign the value zero if vertex too far form center
    Faces2center[Faces2center.>radius_search^2] .= 0
    # % Sum the distance of vertices for each faces
    Faces2center = Faces2center[:,1]+Faces2center[:,2]+Faces2center[:,3]
    # % Create coefficient matrix for neighbourg of center of considered face.
    # % Only faces with non zero distances are valid.
    index_row = find_nonzero_index(Faces2center)

    F_neighbourg[i,1:length(index_row)] = index_row'
    F_neighbourg[i,1+length(index_row)] = i

    if length(index_row)+1 > maximum_neighbour
        maximum_neighbour = length(index_row)+1
    end

end
F_neighbourg[F_neighbourg .== 0] .= NaN
# F_neighbourg[:,maximum_neighbour+1:end] = []  # TODO: checken, was das machen soll


##############################################################
# create folder for plots
##############################################################
if !isdir(folder_plots)
    mkdir(folder_plots)
end
# create folder for particle simulation
if !isdir(folder_particle_simula)
    mkdir(folder_particle_simula)
end


##############################################################
# % Begin - Initialize particles with random initial orientation
##############################################################

r = zeros(num_part,3); # Position of particles
n = zeros(num_part,3); # Orientation of particles

p1 = Vertex[Faces[:,1],:] # Vertex number 1 of each faces
p2 = Vertex[Faces[:,2],:] # Vertex number 2 of each faces
p3 = Vertex[Faces[:,3],:] # Vertex number 3 of each faces
Norm_vect = ones(num_part,3); # Initialisation of normal vector of each faces
num_face = fill(NaN, num_part, length(F_neighbourg[1,:])^2)   # Initialisation of face on which is the particle

# for-loop where particle are randomly positioned and orientated in a 3D
# box, then projected on the closest face. The normal vector of that face
# is then also saved
for i=1:num_part

    # Randomly position and orientate particle number i
    r[i,:] = [(maximum(Vertex[:,1])-minimum(Vertex[:,1]))*rand(1)[1]+minimum(Vertex[:,1]),
        (maximum(Vertex[:,2])-minimum(Vertex[:,2]))*rand(1)[1]+minimum(Vertex[:,2]),
        (maximum(Vertex[:,3])-minimum(Vertex[:,3]))*rand(1)[1]+minimum(Vertex[:,3])
        ]

    #random particle orientation
    n[i,:] = [-1+2*rand(1)[1],-1+2*rand(1)[1],-1+2*rand(1)[1]]

    # Initialize number faces candidates
    face_candi = fill(NaN, 1, length(num_face[i,:]))
    candi_num = 0
    # Vector of particle to faces
    face_coord_temp = Faces_coord - cat(ones(size(Faces_coord[:,:,3]))*
        r[i,1],ones(size(Faces_coord[:,:,3]))*
        r[i,2],ones(size(Faces_coord[:,:,3]))*
        r[i,3], dims=3)

    # Distance of particle to faces
    Dist2faces = sqrt.(sum(face_coord_temp.^2,dims=3)[:,:,1])

    # Faces with closest point to r[i,:] and associated normal vectors
    index_binary = get_index_binary(Dist2faces)
    face_coord_temp = Faces_coord[index_binary,:,:]
    N_temp = N[index_binary,:] |> vec_of_vec_to_array  # transform the vec of vec into array

    p0s = get_particle_position(face_coord_temp, N_temp, i)

    # Check what face in which the projection is
    p1p2 = reshape(face_coord_temp[:,2,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
    p1p3 = reshape(face_coord_temp[:,3,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
    p1p0 = p0s-reshape(face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
    p1p3crossp1p2 = [p1p3[:,2].*p1p2[:,3]-p1p3[:,3].*p1p2[:,2],-p1p3[:,1].*p1p2[:,3]+p1p3[:,3].*p1p2[:,1],p1p3[:,1].*p1p2[:,2]-p1p3[:,2].*p1p2[:,1]] |> vec_of_vec_to_array
    p1p3crossp1p0 = [p1p3[:,2].*p1p0[:,3]-p1p3[:,3].*p1p0[:,2],-p1p3[:,1].*p1p0[:,3]+p1p3[:,3].*p1p0[:,1],p1p3[:,1].*p1p0[:,2]-p1p3[:,2].*p1p0[:,1]] |> vec_of_vec_to_array
    p1p2crossp1p0 = [p1p2[:,2].*p1p0[:,3]-p1p2[:,3].*p1p0[:,2],-p1p2[:,1].*p1p0[:,3]+p1p2[:,3].*p1p0[:,1],p1p2[:,1].*p1p0[:,2]-p1p2[:,2].*p1p0[:,1]] |> vec_of_vec_to_array
    p1p2crossp1p3 = [p1p2[:,2].*p1p3[:,3]-p1p2[:,3].*p1p3[:,2],-p1p2[:,1].*p1p3[:,3]+p1p2[:,3].*p1p3[:,1],p1p2[:,1].*p1p3[:,2]-p1p2[:,2].*p1p3[:,1]] |> vec_of_vec_to_array

    len_p1p3crossp1p2 = length(p1p3crossp1p2[:,1])

    # "index_face_out" are the row index(es) of face_coord_temp in which a
    # particle cannot be projected. 
    # "index_binary(index_face_in)" are the faces number(s) in which the 
    # particle cannot be projected
    index_face_out = (sum(p1p3crossp1p0.*p1p3crossp1p2,dims=2).<0) .|
        (sum(p1p2crossp1p0.*p1p2crossp1p3,dims=2).<0) .|
        ((sqrt.(sum(p1p3crossp1p0.^2,dims=2))+sqrt.(sum(p1p2crossp1p0.^2,dims=2)))./
        sqrt.(sum(p1p2crossp1p3.^2,dims=2)) .> 1) |> Array |> findall 
    index_face_out= first.(Tuple.(index_face_out))

    # % "index_face_in" are the row index(es) of face_coord_temp in which a
    # % particle can be projected. 
    # "index_binary(index_face_in)" are the faces number(s) in which the 
    # particle can be projected
    index_face_in = setdiff(index_face_out, reshape([1:len_p1p3crossp1p2;], :, 1))

    # If the projections are in no face, take the average projection and
    # normal vectors. Save the faces number used
    if isempty(index_face_in) == 1
        Norm_vect[i,:] = mean(N_temp[index_face_out,:],dims=1)
        r[i,:] = mean(p0s[index_face_out,:],dims=1)
        num_face[i,1:length(index_face_out)] = index_binary[index_face_out]'

    # If the projections are in a face, save its number, normal vector and
    # projected point
    else
        index_face_in = index_face_in[1]  # because first particle can be
        # projected in different faces in this initial projection
        Norm_vect[i,:] = N_temp[index_face_in,:]
        num_face[i,1] = index_binary[index_face_in]
        r[i,:] = p0s[index_face_in,:]
    end
end

# % find faces closer to each points and associated normal vector
# %%%%%%%%%%%%

# %Project the orientation of the corresponding faces using normal vectors
n = P_perp(Norm_vect,n)
n = n./(sqrt.(sum(n.^2,dims=2))*ones(1,3)) # And normalise orientation vector

# particle_info[Symbol("t",0)] = [r,n]
#particle_info = setfield(particle_info,num2str(0),[r,n]);

##############################################################
# End - Initialize particles with random initial orientation
##############################################################



# # %%
# # % Initialize graphical output
# if movie_making == 1
#     movie_file = fullfile(UpFolder,"movies","test.mp4")
#     writerObj = VideoWriter(movie_file,"MPEG-4")  # TODO: Open the video writer object
#     writerObj.FrameRate = fps
#     writerObj.Quality = 100
#     open(writerObj)
# end


# %%
# % for loop over each time step. In the first part of the loop, calculate
# % force, displacement and reorientation for all particles. In the second
# % part, for each particle project them on closest face. In the third part,
# % we sometimes plot the position and orientation of the cells
tt=1  # TODO: remove this line
for tt=1:num_step #number of time steps
    r_modified = [];
    # if loop to change forces and velocity after some time because in
    # first time steps, just repulsive force for homogeneisation of
    # particle repartition
    if tt > 500
        v0 = v0_next;
        k = k_next;
    end

    ##############################################################
    # Part 1 of for loop
    ##############################################################
    
    # % Vector between all particles (1->2; 1->3; 1->4;... 641->1; 641->2;
    # % 641->3; 641->4;...641->640...)
    # % Vect is a 3D matrix for vector between particles, e.g. below
    # % Vect(1,2,1) = r(1,1) - r(2,1)
    # % Vect(1,2,2) = r(1,2) - r(2,2)
    # % Vect(1,2,3) = r(1,3) - r(2,3)
    Vect = cat(dims=3,
        r[:,1]*ones(1,num_part)-ones(num_part,1)*r[:,1]',
        r[:,2]*ones(1,num_part)-ones(num_part,1)*r[:,2]',
        r[:,3]*ones(1,num_part)-ones(num_part,1)*r[:,3]'
        )

    # % cross_Nij is the cross product between orientation vectors
    # % cross_Nij(1,2,1) = n(1,2)*n(2,3)-(n(1,3)*n(2,2))
    # % cross_Nij(1,2,2) = -(n(1,1)*n(2,3))+n(1,3)*n(2,1)
    # % cross_Nij(1,2,3) = n(1,1)*n(2,2)-(n(1,2)*n(2,1))
    cross_Nij = cat(dims=3,
        (n[:,2]*ones(1,num_part)).*(ones(num_part,1)*n[:,3]')-
        (n[:,3]*ones(1,num_part)).*(ones(num_part,1)*n[:,2]'),
        -(n[:,1]*ones(1,num_part)).*(ones(num_part,1)*n[:,3]')+
        (n[:,3]*ones(1,num_part)).*(ones(num_part,1)*n[:,1]'),
        (n[:,1]*ones(1,num_part)).*(ones(num_part,1)*n[:,2]')-
        (n[:,2]*ones(1,num_part)).*(ones(num_part,1)*n[:,1]')
        )

    Distmat = sqrt.(sum(Vect.^2,dims=3))[:,:,1]  # distance of each point with the others
    
    Fij_rep = (-k*(2*sigma_n.-Distmat))./(2*sigma_n)
    Fij_adh = (k_adh*(2*sigma_n.-Distmat))./(2*sigma_n-r_adh)

    # No force if particles too far from each other or if particle itself
    Fij_rep[(Distmat .>= 2*sigma_n) .| (Distmat .== 0)] .= 0 
    Fij_adh[(Distmat .< 2*sigma_n) .| (Distmat .> r_adh) .| (Distmat .== 0)] .= 0

    # Fij is the force between particles
    # Fij(1,2,1) = -k(2sigma_n - norm(r(2,:)-r(1,:))) * (r(1,1)-r(2,1)) / norm(r(2,:)-r(1,:))
    # Fij(1,2,2) = -k(2sigma_n - norm(r(2,:)-r(1,:))) * (r(2,2)-r(2,2)) / norm(r(2,:)-r(1,:))
    # Fij(1,2,3) = -k(2sigma_n - norm(r(2,:)-r(1,:))) * (r(2,3)-r(2,3)) / norm(r(2,:)-r(1,:))
    # if 0 < norm(r(2,:)-r(1,:)) < 2*sigma_n
    Fij = Fij_rep .+ Fij_adh
    Fij = cat(dims=3,Fij,Fij,Fij).*(Vect./(cat(dims=3,Distmat,Distmat,Distmat)))
    
    # % Actual force felt by each particle
    F_track = reshape(sum(replace!(Fij, NaN=>0), dims=1),: ,size(Fij,3)) 
    # % Velocity of each particle
    r_dot = P_perp(Norm_vect,v0.*n+mu.*F_track)

    r_dotTest = 0;
    r_prev = r;  #% save current position
    n_prev = n;  #% save current orientation
    r = r+r_dot*dt; # % calculate next position


    ##############################################################
    # Visceck-type n correction adapted from "Phys. Rev. E 74, 061908"
    ##############################################################

    ncross = cat(dims=2,n[:,2].*r_dot[:,3]-n[:,3].*r_dot[:,2],
        -(n[:,1].*r_dot[:,3]-n[:,3].*r_dot[:,1]),
        n[:,1].*r_dot[:,2]-n[:,2].*r_dot[:,1]) ./
        (sqrt.(sum(r_dot.^2,dims=2))*ones(1,3))

    n_cross_correction = (1/tau)*ncross*dt

    new_n = n-cat(dims=2,n[:,2].*n_cross_correction[:,3]-n[:,3].*n_cross_correction[:,2],
        -(n[:,1].*n_cross_correction[:,3]-n[:,3].*n_cross_correction[:,1]),
        n[:,1].*n_cross_correction[:,2]-n[:,2].*n_cross_correction[:,1])

    n = new_n./(sqrt.(sum(new_n.^2,dims=2))*ones(1,3))

    ##############################################################
    # End Visceck-type n correction
    ##############################################################

    Norm_vect = ones(num_part,3);
    for i = 1:length(r[:,1])
        # Search for faces around the particle before displacement in wich the
        # cell could migrate. Only face with at least one vertex within the
        # zone defined by the particle at its center and of radius r_dot*dt are
        # candidates for projection
        radius_search = sqrt.(sum(r_dot[i,:].^2))*dt

        number_faces = replace!(num_face[i,:], NaN=>0)'
        number_faces = Int.(number_faces)
        number_faces = number_faces[number_faces .!= 0]

        if length(number_faces) == 1
            Faces_numbers1 = F_neighbourg[number_faces,:]
            #Faces_numbers1 = Faces_numbers1[isnan(Faces_numbers1) == 0]
            full_number_faces = Faces_numbers1[replace!(Faces_numbers1, NaN=>0)]
            # full_number_faces = Faces_numbers1[isnan(Faces_numbers1) == 0]  # NOTE: I replaced it without checking
        else
            full_number_faces = number_faces'
            num_face_i = 1
            for num_face_i = 1:length(number_faces)
                Faces_numbers1 = F_neighbourg[Int(number_faces[num_face_i]),:]'
                Faces_numbers1 = replace!(Faces_numbers1, NaN=>0)
                Faces_numbers1 = Int.(Faces_numbers1)
                Faces_numbers1 = Faces_numbers1[Faces_numbers1 .!= 0]'  # remove the 0 values
                full_number_faces = cat(dims=2,full_number_faces,Faces_numbers1)
             end
        end
        full_number_faces = unique(full_number_faces)

        # Faces coordinates
        Faces_coord_temp = Faces_coord[full_number_faces,:,:]
        # % Normal vectors of faces
        NV = N[full_number_faces,:] |> vec_of_vec_to_array
        # % Vector of particle to faces
        face_coord_temp = Faces_coord_temp .- cat(dims=3,ones(size(Faces_coord_temp[:,:,3]))*r[i,1],
            ones(size(Faces_coord_temp[:,:,3]))*r[i,2],
            ones(size(Faces_coord_temp[:,:,3]))*r[i,3])

        # % Distance of particle to faces
        Dist2faces = sqrt.(sum(face_coord_temp.^2,dims=3))[:,:,1]

        # TODO: this could all maybe be a function or class
        # % Faces with closest point to r(i,:) and associated normal vectors
        index_binary = get_index_binary(Dist2faces)
        face_coord_temp = Faces_coord_temp[index_binary,:,:]
        N_temp = NV[index_binary,:]

        p0s = get_particle_position(face_coord_temp, N_temp, i)  # get the position of the particle on the face

        # Check what face in which the projection is
        p1p2 = reshape(face_coord_temp[:,2,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
        p1p3 = reshape(face_coord_temp[:,3,:]-face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
        p1p0 = p0s-reshape(face_coord_temp[:,1,:],length(face_coord_temp[:,1,1]),3);
        p1p3crossp1p2 = [p1p3[:,2].*p1p2[:,3]-p1p3[:,3].*p1p2[:,2],-p1p3[:,1].*p1p2[:,3]+p1p3[:,3].*p1p2[:,1],p1p3[:,1].*p1p2[:,2]-p1p3[:,2].*p1p2[:,1]] |> vec_of_vec_to_array
        p1p3crossp1p0 = [p1p3[:,2].*p1p0[:,3]-p1p3[:,3].*p1p0[:,2],-p1p3[:,1].*p1p0[:,3]+p1p3[:,3].*p1p0[:,1],p1p3[:,1].*p1p0[:,2]-p1p3[:,2].*p1p0[:,1]] |> vec_of_vec_to_array
        p1p2crossp1p0 = [p1p2[:,2].*p1p0[:,3]-p1p2[:,3].*p1p0[:,2],-p1p2[:,1].*p1p0[:,3]+p1p2[:,3].*p1p0[:,1],p1p2[:,1].*p1p0[:,2]-p1p2[:,2].*p1p0[:,1]] |> vec_of_vec_to_array
        p1p2crossp1p3 = [p1p2[:,2].*p1p3[:,3]-p1p2[:,3].*p1p3[:,2],-p1p2[:,1].*p1p3[:,3]+p1p2[:,3].*p1p3[:,1],p1p2[:,1].*p1p3[:,2]-p1p2[:,2].*p1p3[:,1]] |> vec_of_vec_to_array

        len_p1p3crossp1p2 = length(p1p3crossp1p2[:,1])

        # "index_face_out" are the row index(es) of face_coord_temp in which a
        # particle cannot be projected. 
        # "index_binary(index_face_in)" are the faces number(s) in which the 
        # particle cannot be projected
        index_face_out = (sum(p1p3crossp1p0.*p1p3crossp1p2,dims=2).<0) .|
            (sum(p1p2crossp1p0.*p1p2crossp1p3,dims=2).<0) .|
            ((sqrt.(sum(p1p3crossp1p0.^2,dims=2))+sqrt.(sum(p1p2crossp1p0.^2,dims=2)))./
            sqrt.(sum(p1p2crossp1p3.^2,dims=2)) .> 1) |> Array |> findall 
        index_face_out= first.(Tuple.(index_face_out))
    
        # % "index_face_in" are the row index(es) of face_coord_temp in which a
        # % particle can be projected. 
        # "index_binary(index_face_in)" are the faces number(s) in which the 
        # particle can be projected
        index_face_in = setdiff(index_face_out, reshape([1:len_p1p3crossp1p2;], :, 1))

        # % If the projections are in no face, take the average projection and
        # % normal vectors. Save the faces number used
        if isempty(index_face_in) == 1
            Norm_vect[i,:] = mean(N_temp[index_face_out,:],dims=1)
            r[i,:] = mean(p0s[index_face_out,:],dims=1)
            num_face[i,1:length(index_face_out)] = full_number_faces[index_binary[index_face_out]]'

        # % If the projections are in a face, save its number, normal vector and
        # % projected point
        else
            if length(index_face_in) > 1
                Norm_vect[i,:] = mean(N_temp[index_face_in,:],dims=1)
                r[i,:] = mean(p0s[index_face_in,:],dims=1)
                num_face[i,1:length(index_face_in)] = full_number_faces[index_binary[index_face_in]]'  # TODO: this is empty
            else
                Norm_vect[i,:] = N_temp[index_face_in,:]
                num_face[i,1] = full_number_faces[index_binary[index_face_in]]   # TODO: this is empty
                r[i,:] = p0s[index_face_in,:]
            end
        end
    end

    # % find faces closer to each points and associated normal vector
    # %%%%%%%%%%%%
    # %Project the orientation of the corresponding faces using normal vectors
    n = P_perp(Norm_vect,n)
    n = n./(sqrt.(sum(n.^2,dims=2))*ones(1,3)) # And normalise orientation vector

    t=(tt-1)*dt;


    ##############################################################
    # Part 3 of for loop%
    ##############################################################

    # %Graphic output
    if rem(tt,plotstep)==0

        # normalize r_dot for visualization
        nr_dot=[];

        for i=1:num_part
            nr_dot=[nr_dot; r_dot[i,:]/norm(r_dot[i,:])];
        end
        
        nr_dot_cross=[];
        for i=1:num_part
            cross_Nrdot=cross(n[i,:],Norm_vect[i,:])
            nr_dot_cross=[nr_dot_cross; cross_Nrdot./norm(cross_Nrdot)]   # TODO: maybe check this again
        end

        # %evaluate number of neighbourgs within 2.4 sigma cut off
        num_partic = ones(size(Distmat));
        num_partic[(Distmat .== 0) .| (Distmat .> 2.4*sigma_n)] .= 0;
        # %list of nearest neighbour to each particle
        number_neighbours=sum(num_partic,dims=2);

        # % specify colorvalue "c" according to neighbours for particles to plot
        N_color=[];
        for i=1:num_part
            append!(N_color, Int.(number_neighbours[i,:]))
        end

        c=N_color;


        ##############################################################
        # Definition of the order parameter:
        ##############################################################

        # Define a vector normal to position vector and velocity vector
        v_tp=[r[:,2].*r_dot[:,3]-r[:,3].*r_dot[:,2],-r[:,1].*r_dot[:,3]+r[:,3].*r_dot[:,1],r[:,1].*r_dot[:,2]-r[:,2].*r_dot[:,1]];
        v_tp = v_tp |> vec_of_vec_to_array |> transpose
        # Normalize v_tp
        v_norm=v_tp./(sqrt.(sum(v_tp.^2,dims=2))*ones(1,3));
        # Sum v_tp vectors and devide by number of particle to obtain order
        # parameter of collective motion for spheroids
        v_order[tt/plotstep]=(1/num_part)*norm(sum(v_norm,dims=1));

        # #Calculate angles for equirectangular map projection
        # phi1 = asin(r[:,3]/sqrt.(sum(r.^2,dims=2)));   # elevation angle
        # thet1 = atan2(r[:,2],r[:,1]); # azimuth


        #custom colormap for the visualization of nearest neighbours (purple, red, blue, green, yellow)
        # cmap=[0.6, 0.0, 1.0; 1.0, 0.0, 0.0; 0.0, 0.0, 0.8; 0.0, 1.0, 0.0; 1.0, 1.0, 0.0];

        #set figure settings
#         set(hFig, 'Position', [50 10 1280 720], 'color', [1 1 1],'Renderer','OpenGL')
# # %         set(hFig, 'Position', [50 50 1280 720], 'color', [1 1 1])

#          subplot(3,5,[1 13]); #plot the sphere and particles in 3D
#         #--------------------------------------------------------------------------
#         # % Plot the surface
#         [xs,ys,zs]=ellipsoid(0,0,0,maximum(Vertex[:,1]),maximum(Vertex[:,2]),maximum(Vertex[:,3]),30);
#         ss=1;
#         h1=surf(xs,ys,zs);
#         set (h1,'EdgeColor',[0.75,0.75,0.75],'FaceColor',[0.95,0.95,0.75],'MeshStyle','row');
#                 alpha(0.7);
#            camzoom(1.6);
        

#         # % Plot the particles    msize = Observable(initms) # this is the variable that will change
            # first frame, initial plot

#         # %Definition of contourplot to show particle density
#         F = TriScatteredInterp(thet1,phi1,number_neighbours);
#         t_thet1= -pi:2*pi/100:pi;
# # %         t_phi1= -pi/2:pi/100:pi/2;
#         t_phi1= -1.45:pi/100:1.45;
#         [q_thet1,q_phi1]=meshgrid(t_thet1,t_phi1);
#         q_neighbour=F(q_thet1,q_phi1);
#         q_neighbour(q_neighbour==0)=nan;

#         # %squared norm of each element in r_dot
#         ner_dot=[];
#         for i=1:num_part
#             ner_dot=[ner_dot; norm(r_dot(i,:))^2];
#         end

#         # %total energy in the system
#         if 0.5*sum(ner_dot)<200
#             en_v(tt/plotstep)=0.5.*sum(ner_dot);
#         end

#         x_en_v=0:(length(en_v)-1);
#         x_v_order=0:(length(v_order)-1);
    end
end

