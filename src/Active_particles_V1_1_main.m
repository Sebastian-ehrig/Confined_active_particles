clear all;
close all;
tic

addpath('./functions')
UpFolder = fileparts(pwd);

% Import an STL mesh, returning a PATCH-compatible face-vertex structure
namestructure = 'ellipsoid_x4_y2.5_z1_coarse_mesh_refSph6_sigma_0.416_N829_Area452.389';

%namestructure = 'prolate_spheroid_aspr8_coarse_mesh_area621.0128area_refSph7';
[Faces,Vertex,N] = stlread(fullfile(UpFolder,'meshes',[namestructure,'.stl']));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When we plot surfaces, we are interested only in some part of the
% geometry, so we select the faces which have at least one vertex in the
% area of interest
V_plot=ones(length(Faces(:,1)),1); % Will store the num of the faces in the ROI
face_num_plot=[];
Faces_coord = cat(3,[Vertex(Faces(:,1),1),Vertex(Faces(:,2),1),...
    Vertex(Faces(:,3),1)],[Vertex(Faces(:,1),2),Vertex(Faces(:,2),2),...
    Vertex(Faces(:,3),2)],[Vertex(Faces(:,1),3),Vertex(Faces(:,2),3),...
    Vertex(Faces(:,3),3)]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Define Model Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visible_off = 0; % 0 to display figure while plotting them, 1 to just save 
                 % them without being displayed first
                 
movie_making = 1; % 1 for movies, 0 for saving images

num_part = 828; % number particles
v0 = 0.1; % velocity of particles
v0_next = 0.1;% if velocity is to be changed after a certain number of timesteps

num_step =  3000; %Number of timesteps

Area=452.389; % total surface area (should be extracted from mesh-name)

sigma_n = 5/12;% particle radius

phi_pack=(num_part*pi*(sigma_n^2))/Area; % packing-fraction (~density)

% Define parameters for the force calculations:
% ---------------------------------------------
r_adh = 1; % Cut off radius for adhesive forces
k = 10.; % Elastic constant for repulsive forces
k_next = 10.; % value of elastic constant after a certain number of timesteps
k_adh = 0.75; % Adhesive forces

mu = 1.; % arbitrary mass of particle
tau = 1; %time of relaxation of collision between 2 particles

dt = 0.01*tau; % Step size
plotstep = 0.1/dt; % Number of calculation between plot

% a = 0.5; % arbitrary value
b = v0_next/(2*sigma_n*mu*k_next); % for calculating coupling constant
J = b/tau; % Coupling constant of orientation vectors
noise=0.002/tau;  %value of the noise amplitude
hi=noise/(2*sqrt(dt)); % Amplitude of the noise on orientation

% Parameters for making movies:
%------------------------------
rho = 6; % arbitrary value for scaling the vector on plotting. To automatize, maybe from file name value
scale=0.1*rho; % Scale the vector for making movies
fps = 25; % Frame per second of movie



v_order = zeros(num_step/plotstep,1);
en_v = zeros(num_step/plotstep,1);
tot_F = zeros(num_step/plotstep,1);

particle_info = struct;


name_data = [namestructure,'_N_',num2str(num_part),'_rho_',num2str(rho)...
        ,'_F_rep_',num2str(k),'_F_adh_',num2str(k_adh),'_tau_',num2str(tau),'_b_',num2str(b),...
        '_s',num2str(v0),'_t',num2str(num_step),'_fps',num2str(25)];
    
% Define folder structure and pre-process meshes:
% -----------------------------------------------
mesh_struct = fullfile(UpFolder,'meshes',[namestructure,'_proj.mat']);
folder_plots = fullfile(UpFolder,'images',namestructure);
folder_particle_simula = fullfile(UpFolder,'simulation_structure');
if exist(mesh_struct,'file') == 2
    load(mesh_struct);
    F_neighbourg = F.mesh;
else
    F_neighbourg = Mesh_neighbour_projV1_1(namestructure,1.1*v0);
end
if exist(folder_plots,'dir') ~= 7
     mkdir(folder_plots)
end
if exist(folder_particle_simula,'dir') ~= 7
     mkdir(folder_particle_simula)
end

%%
%%%%%%%%%%%%%%%%%
% Define perpendicular projection functions, adapted from "PhysRevE 91 022306"
%%%%%%%%%%%%%%%%%

P_perp = @(a,b) (b-(sum(b.*a,2)./(sqrt(sum(a.^2,2)).^2)*ones(1,3)).*a);
% % P_perp does a normal projection of the vector b on the plane normal to a
% % , the output c is a 2d array, num_part x 3
P_plan = @(a,b,a1) ((sum(b.*a,2)./(sqrt(sum(a.^2,2)))*ones(1,3)).*...
    [a(:,2).*a1(:,3)-a(:,3).*a1(:,2),-a(:,1).*a1(:,3)+a(:,3).*a1(:,1),...
    a(:,1).*a1(:,2)-a(:,2).*a1(:,1)]);
% % P_plan does a projection of the vector b on vector normal to a1 in
% %the plane normal to a, the output c is a 2d array, num_part x 3

%%
%%%%%%%%%%%%%%%%%
% Begin - Initialize particles with random initial orientation
%%%%%%%%%%%%%%%%%%

r = zeros(num_part,3); % Position of particles
n = zeros(num_part,3); % Orientation of particles

p1 = Vertex(Faces(:,1),:); % Vertex number 1 of each faces
p2 = Vertex(Faces(:,2),:); % Vertex number 2 of each faces
p3 = Vertex(Faces(:,3),:); % Vertex number 3 of each faces
Norm_vect = ones(num_part,3); % Initialisation of normal vector of each faces
num_face = NaN(num_part,length(F_neighbourg(1,:))^2); % Initialisation of face on which is the particle
% for-loop where particle are randomly positioned and orientated in a 3D
% box, then projected on the closest face. The normal vector of that face
% is then also saved
for i=1:num_part
    % Randomly position and orientate particle number i
    r(i,:) = [(max(Vertex(:,1))-min(Vertex(:,1)))*rand(1)+min(Vertex(:,1)),...
        (max(Vertex(:,2))-min(Vertex(:,2)))*rand(1)+min(Vertex(:,2)),...
        (max(Vertex(:,3))-min(Vertex(:,3)))*rand(1)+min(Vertex(:,3))];
    
%     r(i,:)=[((1+i)*sigma_n/i),((1+i)*sigma_n/i),((1+i)*sigma_n/i)];
%     n(i,:) = [1,1,1];

    %random particle orientation
    n(i,:) = [-1+2.*rand(1),-1+2.*rand(1),-1+2.*rand(1)];
%    

   
    % Initialize number faces candidates
    face_candi = NaN(1,length(num_face(i,:)));
    candi_num = 0;
    % Vector of particle to faces
    face_coord_temp = Faces_coord - cat(3,ones(size(Faces_coord(:,:,3)))*...
        r(i,1),ones(size(Faces_coord(:,:,3)))*...
        r(i,2),ones(size(Faces_coord(:,:,3)))*...
        r(i,3));
    % Distance of particle to faces
    Dist2faces = sqrt(sum(face_coord_temp.^2,3));
    % Faces with closest point to r(i,:) and associated normal vectors
    index_binary = find(sum(Dist2faces == min(min(Dist2faces)),2) > 0);
    face_coord_temp = Faces_coord(index_binary,:,:);
    N_temp = N(index_binary,:);
    particle = cat(2,ones(size(face_coord_temp(:,1,3)))*...
        r(i,1),ones(size(face_coord_temp(:,1,3)))*...
        r(i,2),ones(size(face_coord_temp(:,1,3)))*...
        r(i,3));
    p0s = particle - cat(2,sum((particle-reshape(face_coord_temp(:,1,:),...
        length(face_coord_temp(:,1,1)),3)).*N_temp,2),...
        sum((particle-reshape(face_coord_temp(:,1,:),...
        length(face_coord_temp(:,1,1)),3)).*N_temp,2),...
        sum((particle-reshape(face_coord_temp(:,1,:),...
        length(face_coord_temp(:,1,1)),3)).*N_temp,2)...
        ).*N_temp;
    % Check what face in which the projection is
    p1p2 = reshape(face_coord_temp(:,2,:)-face_coord_temp(:,1,:),...
        length(face_coord_temp(:,1,1)),3);
    p1p3 = reshape(face_coord_temp(:,3,:)-face_coord_temp(:,1,:),...
        length(face_coord_temp(:,1,1)),3);
    p1p3crossp1p2 = [p1p3(:,2).*p1p2(:,3)-p1p3(:,3).*p1p2(:,2),-p1p3(:,1).*...
        p1p2(:,3)+p1p3(:,3).*p1p2(:,1),p1p3(:,1).*p1p2(:,2)-p1p3(:,2).*p1p2(:,1)];
    p1p0 = p0s-reshape(face_coord_temp(:,1,:),length(face_coord_temp(:,1,1)),...
        3);
    p1p3crossp1p0 = [p1p3(:,2).*p1p0(:,3)-p1p3(:,3).*p1p0(:,2),...
        -p1p3(:,1).*p1p0(:,3)+p1p3(:,3).*p1p0(:,1),...
        p1p3(:,1).*p1p0(:,2)-p1p3(:,2).*p1p0(:,1)];
    
    p1p2crossp1p0 = [p1p2(:,2).*p1p0(:,3)-p1p2(:,3).*p1p0(:,2),...
        -p1p2(:,1).*p1p0(:,3)+p1p2(:,3).*p1p0(:,1),...
        p1p2(:,1).*p1p0(:,2)-p1p2(:,2).*p1p0(:,1)];
    p1p2crossp1p3 = [p1p2(:,2).*p1p3(:,3)-p1p2(:,3).*p1p3(:,2),...
        -p1p2(:,1).*p1p3(:,3)+p1p2(:,3).*p1p3(:,1),...
        p1p2(:,1).*p1p3(:,2)-p1p2(:,2).*p1p3(:,1)];
    
    % "index_face_in" are the row index(es) of face_coord_temp in which a
    % particle can be projected. 
    % "index_binary(index_face_in)" are the faces number(s) in which the 
    % particle can be projected
    index_face_in = find((sum(p1p3crossp1p0.*p1p3crossp1p2,2)>=0) & (sum(p1p2crossp1p0.*p1p2crossp1p3,2)>=0)...
        & ((sqrt(sum(p1p3crossp1p0.^2,2))+sqrt(sum(p1p2crossp1p0.^2,2)))./...
        sqrt(sum(p1p2crossp1p3.^2,2)) <= 1));
    % "index_face_out" are the row index(es) of face_coord_temp in which a
    % particle cannot be projected. 
    % "index_binary(index_face_in)" are the faces number(s) in which the 
    % particle cannot be projected
    index_face_out = find((sum(p1p3crossp1p0.*p1p3crossp1p2,2)<0) | (sum(p1p2crossp1p0.*p1p2crossp1p3,2)<0)...
        | ((sqrt(sum(p1p3crossp1p0.^2,2))+sqrt(sum(p1p2crossp1p0.^2,2)))./...
        sqrt(sum(p1p2crossp1p3.^2,2)) > 1));
    
    % If the projections are in no face, take the average projection and
    % normal vectors. Save the faces number used
    if isempty(index_face_in) == 1
        Norm_vect(i,:) = mean(N_temp(index_face_out,:),1); 
        r(i,:) = mean(p0s(index_face_out,:),1);
        num_face(i,1:length(index_face_out)) = index_binary(index_face_out)';
    % If the projections are in a face, save its number, normal vector and
    % projected point
    else
        index_face_in = index_face_in(1); % because first particle can be
        % projected in different faces in this initial projection
        Norm_vect(i,:) = N_temp(index_face_in,:);
        num_face(i,1) = index_binary(index_face_in);
        r(i,:) = p0s(index_face_in,:);
    end
end
% find faces closer to each points and associated normal vector
%%%%%%%%%%%%

%Project the orientation of the corresponding faces using normal vectors
n = P_perp(Norm_vect,n);
n = n./(sqrt(sum(n.^2,2))*ones(1,3)); % And normalise orientation vector
particle_info.(['t',num2str(0)]) = [r,n];
%particle_info = setfield(particle_info,num2str(0),[r,n]);

%%%%%%%%%%%%%%%%%
% End - Initialize particles with random initial orientation
%%%%%%%%%%%%%%%%%%

%%
% Initialize graphical output
if movie_making == 1
    movie_file = fullfile(UpFolder,'movies',[name_data,'.mp4']);
    writerObj = VideoWriter(movie_file,'MPEG-4'); % Open the video writer object
    writerObj.FrameRate = fps;
    writerObj.Quality = 100;
    open(writerObj);
end

%%
% for loop over each time step. In the first part of the loop, calculate
% force, displacement and reorientation for all particles. In the second
% part, for each particle project them on closest face. In the third part,
% we sometimes plot the position and orientation of the cells

for tt=1:num_step %number of time steps
    r_modified = [];
    % if loop to change forces and velocity after some time because in
    % first time steps, just repulsive force for homogeneisation of
    % particle repartition
    if tt > 500
        v0 = v0_next;
        k = k_next;
%     elseif tt==500 && v0<0.1
%         n = -1+2*rand(num_part,3);
%         n = P_perp(Norm_vect,n);
%         n = n./(sqrt(sum(n.^2,2))*ones(1,3)); % And normalise orientation vector
    end
    
    %%%%%%%%%%%%%%
    % Part 1 of for loop
    %%%%%%%%%%%%%%
    
    % Vector between all particles (1->2; 1->3; 1->4;... 641->1; 641->2;
    % 641->3; 641->4;...641->640...)
    Vect = cat(3,...
        r(:,1)*ones(1,num_part)-ones(num_part,1)*r(:,1)',...
        r(:,2)*ones(1,num_part)-ones(num_part,1)*r(:,2)',...
        r(:,3)*ones(1,num_part)-ones(num_part,1)*r(:,3)'...
        );
    % Vect is a 3D matrix for vector between particles, e.g. below
    % Vect(1,2,1) = r(1,1) - r(2,1)
    % Vect(1,2,2) = r(1,2) - r(2,2)
    % Vect(1,2,3) = r(1,3) - r(2,3)
    cross_Nij = cat(3,...
        (n(:,2)*ones(1,num_part)).*(ones(num_part,1)*n(:,3)')-...
        (n(:,3)*ones(1,num_part)).*(ones(num_part,1)*n(:,2)'),...
        -(n(:,1)*ones(1,num_part)).*(ones(num_part,1)*n(:,3)')+...
        (n(:,3)*ones(1,num_part)).*(ones(num_part,1)*n(:,1)'),...
        (n(:,1)*ones(1,num_part)).*(ones(num_part,1)*n(:,2)')-...
        (n(:,2)*ones(1,num_part)).*(ones(num_part,1)*n(:,1)'));
    % cross_Nij is the cross product between orientation vectors
    % cross_Nij(1,2,1) = n(1,2)*n(2,3)-(n(1,3)*n(2,2))
    % cross_Nij(1,2,2) = -(n(1,1)*n(2,3))+n(1,3)*n(2,1)
    % cross_Nij(1,2,3) = n(1,1)*n(2,2)-(n(1,2)*n(2,1))
    
    Distmat=sqrt(sum(Vect.^2,3)); %distance of each point with the others
    
    Fij_rep = (-k*(2*sigma_n-Distmat))./(2*sigma_n);
    Fij_rep((Distmat >= 2*sigma_n) | (Distmat == 0))= 0; % No force if
    % particles too far from each other or if particle itself
    
    Fij_adh = (k_adh*(2*sigma_n-Distmat))./(2*sigma_n-r_adh);
    Fij_adh((Distmat < 2*sigma_n) | (Distmat > r_adh) | (Distmat == 0))= 0; % No force if
    % particles too far from each other or if particle itself
    
    Fij = Fij_rep+Fij_adh;
    Fij = cat(3,Fij,Fij,Fij).*(Vect./(cat(3,Distmat,Distmat,Distmat)));
    % Fij is the force between particles
    % Fij(1,2,1) = -k(2sigma_n - norm(r(2,:)-r(1,:))) * (r(1,1)-r(2,1)) / norm(r(2,:)-r(1,:))
    % Fij(1,2,2) = -k(2sigma_n - norm(r(2,:)-r(1,:))) * (r(2,2)-r(2,2)) / norm(r(2,:)-r(1,:))
    % Fij(1,2,3) = -k(2sigma_n - norm(r(2,:)-r(1,:))) * (r(2,3)-r(2,3)) / norm(r(2,:)-r(1,:))
    % if 0 < norm(r(2,:)-r(1,:)) < 2*sigma_n
    
    % Actual force felt by each particle
    F_track = reshape(nansum(Fij,1),[],size(Fij,3));
    % Velocity of each particle
    r_dot = P_perp(Norm_vect,v0*n+mu*reshape(nansum(Fij,1),[],size(Fij,3)));
    
    r_dotTest = 0;
    r_prev = r; % save current position
    n_prev = n; % save current orientation
    r = r+r_dot*dt; % calculate next position
    
    %%%%%%%%%%%%%%%%%%%%%
    % Visceck-type n correction adapted from "Phys. Rev. E 74, 061908"
    %%%%%%%%%%%%%%%%%%%%%%
    ncross = cat(2,n(:,2).*r_dot(:,3)-n(:,3).*r_dot(:,2),...
        -(n(:,1).*r_dot(:,3)-n(:,3).*r_dot(:,1)),...
        n(:,1).*r_dot(:,2)-n(:,2).*r_dot(:,1))./(sqrt(sum(r_dot.^2,2))*ones(1,3));
    n_cross_correction = (1/tau)*ncross*dt;
    new_n = n-cat(2,n(:,2).*n_cross_correction(:,3)-n(:,3).*n_cross_correction(:,2),...
        -(n(:,1).*n_cross_correction(:,3)-n(:,3).*n_cross_correction(:,1)),...
        n(:,1).*n_cross_correction(:,2)-n(:,2).*n_cross_correction(:,1)); %n+cross(n_cross_correction,n)
    n = new_n./(sqrt(sum(new_n.^2,2))*ones(1,3));
    %%%%%%%%%%%%%%%%%%%%%
    % End Visceck-type n correction
    %%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%
    % Part 2 of for loop%
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%
    
    Norm_vect = ones(num_part,3);
    for i = 1:length(r(:,1))
        % Search for faces around the particle before displacement in wich the
        % cell could migrate. Only face with at least one vertex within the
        % zone defined by the particle at its center and of radius r_dot*dt are
        % candidates for projection
        radius_search = sqrt(sum(r_dot(i,:).^2))*dt;
        number_faces = num_face(i,isnan(num_face(i,:)) == 0);
        if length(number_faces) == 1
            Faces_numbers1 = F_neighbourg(number_faces,:);
            %Faces_numbers1 = Faces_numbers1(isnan(Faces_numbers1) == 0);
            full_number_faces = Faces_numbers1(isnan(Faces_numbers1) == 0);
        else
            full_number_faces = number_faces;
            for num_face_i = 1:length(number_faces)
                Faces_numbers1 = F_neighbourg(number_faces(num_face_i),:);
                Faces_numbers1 = Faces_numbers1(isnan(Faces_numbers1) == 0);
                full_number_faces = cat(2,full_number_faces,Faces_numbers1);
             end
        end
        full_number_faces = unique(full_number_faces);
        Faces_coord_temp = Faces_coord(full_number_faces,:,:);
        NV = N(full_number_faces,:);
        % Vector of particle to faces
        face_coord_temp = Faces_coord_temp - cat(3,ones(size(Faces_coord_temp(:,:,3)))*...
            r(i,1),ones(size(Faces_coord_temp(:,:,3)))*...
            r(i,2),ones(size(Faces_coord_temp(:,:,3)))*...
            r(i,3));
        % Distance of particle to faces
        Dist2faces = sqrt(sum(face_coord_temp.^2,3));
        % Faces with closest point to r(i,:) and associated normal vectors
        index_binary = find(sum(Dist2faces == min(min(Dist2faces)),2) > 0);
        face_coord_temp = Faces_coord_temp(index_binary,:,:);
        N_temp = NV(index_binary,:);
        
        particle = cat(2,ones(size(face_coord_temp(:,1,3)))*...
            r(i,1),ones(size(face_coord_temp(:,1,3)))*...
            r(i,2),ones(size(face_coord_temp(:,1,3)))*...
            r(i,3));
        p0s = particle - cat(2,sum((particle-reshape(face_coord_temp(:,1,:),...
            length(face_coord_temp(:,1,1)),3)).*N_temp,2),...
            sum((particle-reshape(face_coord_temp(:,1,:),...
            length(face_coord_temp(:,1,1)),3)).*N_temp,2),...
            sum((particle-reshape(face_coord_temp(:,1,:),...
            length(face_coord_temp(:,1,1)),3)).*N_temp,2)...
            ).*N_temp;
        % Check what face in which the projection is
        p1p2 = reshape(face_coord_temp(:,2,:)-face_coord_temp(:,1,:),...
            length(face_coord_temp(:,1,1)),3);
        p1p3 = reshape(face_coord_temp(:,3,:)-face_coord_temp(:,1,:),...
            length(face_coord_temp(:,1,1)),3);
        p1p3crossp1p2 = [p1p3(:,2).*p1p2(:,3)-p1p3(:,3).*p1p2(:,2),-p1p3(:,1).*...
            p1p2(:,3)+p1p3(:,3).*p1p2(:,1),p1p3(:,1).*p1p2(:,2)-p1p3(:,2).*p1p2(:,1)];
        p1p0 = p0s-reshape(face_coord_temp(:,1,:),length(face_coord_temp(:,1,1)),...
            3);
        p1p3crossp1p0 = [p1p3(:,2).*p1p0(:,3)-p1p3(:,3).*p1p0(:,2),...
            -p1p3(:,1).*p1p0(:,3)+p1p3(:,3).*p1p0(:,1),...
            p1p3(:,1).*p1p0(:,2)-p1p3(:,2).*p1p0(:,1)];
        
        p1p2crossp1p0 = [p1p2(:,2).*p1p0(:,3)-p1p2(:,3).*p1p0(:,2),...
            -p1p2(:,1).*p1p0(:,3)+p1p2(:,3).*p1p0(:,1),...
            p1p2(:,1).*p1p0(:,2)-p1p2(:,2).*p1p0(:,1)];
        p1p2crossp1p3 = [p1p2(:,2).*p1p3(:,3)-p1p2(:,3).*p1p3(:,2),...
            -p1p2(:,1).*p1p3(:,3)+p1p2(:,3).*p1p3(:,1),...
            p1p2(:,1).*p1p3(:,2)-p1p2(:,2).*p1p3(:,1)];
        
        % "index_face_in" are the row index(es) of face_coord_temp in which a
        % particle can be projected.
        % "index_binary(index_face_in)" are the faces number(s) in which the
        % particle can be projected
        index_face_in = find((sum(p1p3crossp1p0.*p1p3crossp1p2,2)>=0) & (sum(p1p2crossp1p0.*p1p2crossp1p3,2)>=0)...
            & ((sqrt(sum(p1p3crossp1p0.^2,2))+sqrt(sum(p1p2crossp1p0.^2,2)))./...
            sqrt(sum(p1p2crossp1p3.^2,2)) <= 1));
        % "index_face_out" are the row index(es) of face_coord_temp in which a
        % particle cannot be projected.
        % "index_binary(index_face_in)" are the faces number(s) in which the
        % particle cannot be projected
        index_face_out = find((sum(p1p3crossp1p0.*p1p3crossp1p2,2)<0) | (sum(p1p2crossp1p0.*p1p2crossp1p3,2)<0)...
            | ((sqrt(sum(p1p3crossp1p0.^2,2))+sqrt(sum(p1p2crossp1p0.^2,2)))./...
            sqrt(sum(p1p2crossp1p3.^2,2)) > 1));
        
        % If the projections are in no face, take the average projection and
        % normal vectors. Save the faces number used
        if isempty(index_face_in) == 1
            Norm_vect(i,:) = mean(N_temp(index_face_out,:),1);
            r(i,:) = mean(p0s(index_face_out,:),1);
            num_face(i,1:length(index_face_out)) = full_number_faces(index_binary(index_face_out))';
            % If the projections are in a face, save its number, normal vector and
            % projected point
        else
            if length(index_face_in) > 1
                Norm_vect(i,:) = mean(N_temp(index_face_in,:),1);
                r(i,:) = mean(p0s(index_face_in,:),1);
                num_face(i,1:length(index_face_in)) = full_number_faces(index_binary(index_face_in))';
            else
                Norm_vect(i,:) = N_temp(index_face_in,:);
                num_face(i,1) = full_number_faces(index_binary(index_face_in));
                r(i,:) = p0s(index_face_in,:);
            end
        end
    end
    % find faces closer to each points and associated normal vector
    %%%%%%%%%%%%
    n = P_perp(Norm_vect,n);
    n = n./(sqrt(sum(n.^2,2))*ones(1,3));
    t=(tt-1)*dt;
    
    %particle_info.num2str(tt) = [r,n];
    particle_info.(['t',num2str(tt)]) = [r,n];
 %%   
    %%%%%%%%%%%%%%%%%%%%%
    % Part 3 of for loop%
    %%%%%%%%%%%%%%%%%%%%%  
    
    %Graphic output
    
    if rem(tt,plotstep)==0
        
%       normalize r_dot for visualization
        nr_dot=[];
        for i=1:num_part
        nr_dot=[nr_dot; r_dot(i,:)/norm(r_dot(i,:))];
        end
        
        nr_dot_cross=[];
        for i=1:num_part
            
%         ncross2=[Norm_vect(i,2).*r_dot(i,3)-Norm_vect(i,3).*r_dot(i,2),...
%         -(Norm_vect(i,1).*r_dot(i,3)-Norm_vect(i,3).*r_dot(i,1)),...
%         Norm_vect(i,1).*r_dot(i,2)-Norm_vect(i,2).*r_dot(i,1)];
        cross_Nrdot=cross(n(i,:),Norm_vect(i,:));
        nr_dot_cross=[nr_dot_cross; cross_Nrdot./norm(cross_Nrdot)];
        end
        
        

        %evaluate number of neighbourgs within 2.4 sigma cut off
        num_partic = ones(size(Distmat));
        num_partic((Distmat == 0)|(Distmat > 2.4*sigma_n)) = 0;
        %list of nearest neighbour to each particle
        number_neighbours=sum(num_partic,2);
        
        % specify colorvalue "c" according to forces for particles to plot
%         F_color=[]; 
%         for i=1:num_part
%         F_color=[F_color norm(sum(F_track(i,:)))];
%         end
%         c=F_color;
        % specify colorvalue "c" according to neighbours for particles to plot
        N_color=[]; 
        for i=1:num_part
        %F_color=[F_color norm(sum(F_track(i,:)))];
        N_color=[N_color number_neighbours(i,:)];
        end
        c=N_color;
        % Definition of the order parameter:
        %-----------------------------------
        % Define a vector normal to position vector and velocity vector
        v_tp=[r(:,2).*r_dot(:,3)-r(:,3).*r_dot(:,2),-r(:,1).*r_dot(:,3)+r(:,3).*r_dot(:,1),...
            r(:,1).*r_dot(:,2)-r(:,2).*r_dot(:,1)];
        % Normalize v_tp
        v_norm=v_tp./(sqrt(sum(v_tp.^2,2))*ones(1,3));
        % Sum v_tp vectors and devide by number of particle to obtain order
        % parameter of collective motion for spheroids
        v_order(tt/plotstep)=(1/num_part)*norm(sum(v_norm,1));
        
        %Calculate angles for equirectangular map projection
        phi1 = asin(r(:,3)./sqrt(sum(r.^2,2)));   % elevation angle
        thet1 = atan2(r(:,2),r(:,1)); % azimuth
        
%         phi1_dummy_points=-pi/2:pi/2/10:pi/2;
%         thet1_dummy_points=-pi:pi/10:pi;
%         phi1= [phi1; phi1_dummy_points'];
%         thet1= [thet1; thet1_dummy_points'];
        %gall-peters projection
%         phi1 = 2*sin(phi1);   

        %for cylinder projection
%         phi1=r(:,3); %eleavation
%         thet1 = atan2(r(:,2),r(:,1)); % azimuth

        

        %Draw figure
        hFig = figure(1); 
        if visible_off == 1
            set(hFig, 'Visible', 'off');
        end
        %custom colormap for the visualization of nearest neighbours (purple, red, blue, green, yellow)
        cmap=[0.6, 0.0, 1.0; 1.0, 0.0, 0.0; 0.0, 0.0, 0.8; 0.0, 1.0, 0.0; 1.0, 1.0, 0.0];

        %set figure settings
        set(hFig, 'Position', [50 10 1280 720], 'color', [1 1 1],'Renderer','OpenGL')
%         set(hFig, 'Position', [50 50 1280 720], 'color', [1 1 1])

         subplot(3,5,[1 13]); %plot the sphere and particles in 3D
        %--------------------------------------------------------------------------
       % Plot the surface
        [xs,ys,zs]=ellipsoid(0,0,0,max(Vertex(:,1)),max(Vertex(:,2)),max(Vertex(:,3)),30);
        ss=1;
        h1=surf(xs,ys,zs);
        grid off
        set (h1,'EdgeColor',[0.75,0.75,0.75],'FaceColor',[0.95,0.95,0.75],'MeshStyle','row');
                alpha(0.7);
           camzoom(1.6);

        
%         faceColor  = [0.95,0.95,0.75];
%         edgeColor = [0.75,0.75,0.75];
% %         hold on     
%         P=patch('faces', Faces, 'vertices', Vertex);
%         set(P, 'Facecolor', faceColor,'Facealpha',0.7)
%         set(P, 'Edgecolor', edgeColor, 'Edgealpha',0.8)
        
        
%%%%%%%%%%%%%%%
%         tri = delaunay(Vertex(:,1),Vertex(:,2),Vertex(:,3));
%         h1=trimesh(tri,Vertex(:,1),Vertex(:,2),Vertex(:,3));
%         set (h1,'EdgeColor',[0.75,0.75,0.75],'FaceColor',[0.95,0.95,0.75],'facealpha',0.5);
%%%%%%%%%%%%%%%%%
        hold on
        
        
        %particle positions color coded according to forces:
%         plot3k({r(:,1) r(:,2) r(:,3)},'ColorData',c(:),'ColorRange',[3.5 8.5],'Marker',{'o',6},'ColorBar',false);
%         plot3k({r(:,1) r(:,2) r(:,3)},'ColorData',c(:),'Marker',{'o',6},'ColorBar',false);
        
        scatter3(r(:,1),r(:,2),r(:,3),'filled','CData',c(:))
%         camzoom(1.2);

%            camzoom(1.6);
           
        %position color bar
        pos=get(gca,'pos');
        set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*0.95]);
        pos=get(gca,'pos');
        hc=colorbar('location','northoutside','position',[0.5*pos(1) 0.5*pos(2)+pos(4) 0.5*pos(3) 0.03]);
        set(hc,'xaxisloc','top');
        %axis range
        xlabel(hc, 'No. of nearest neighbour')
%         xlabel(hc, 'force magnitude F/F_{0}')
        grid on
        colormap(cmap);

        caxis([3.5 8.5]);
%         caxis([4.5 7.5]);

        
%         caxis([0 1]);
        hold on
%         
%         r=gpuArray(r);
%         n=gpuArray(n);
        %Draw position of Dummy particle
%         plot3(r(1,1),r(1,2),r(1,3),'o','MarkerSize', 6,'MarkerEdgeColor','k','MarkerFaceColor','r');
%         r=gpuArray(r);
        
%         hold on
%---------------------------------------------                
%         grid on
%         axis([0.67*min(Vertex(:,3))  0.67*max(Vertex(:,3))  0.67*min(Vertex(:,3))  0.67*max(Vertex(:,3))  min(Vertex(:,3))  max(Vertex(:,3))]);
%         axis([0.67*min(Vertex(:,2))-1  0.67*max(Vertex(:,2))+1  0.67*min(Vertex(:,2))-1  0.67*max(Vertex(:,2))+1  min(Vertex(:,2))-1  max(Vertex(:,2))+1]);
%         axis([0.67*min(Vertex(:,3))  0.67*max(Vertex(:,3))  0.67*min(Vertex(:,3))  0.67*max(Vertex(:,3))  min(Vertex(:,3))+0.2*max(Vertex(:,3))  max(Vertex(:,3))-0.2*max(Vertex(:,3))]);
          %axis([0.67*min(Vertex(:,3))  0.67*max(Vertex(:,3))  0.67*min(Vertex(:,3))  0.67*max(Vertex(:,3))  min(Vertex(:,3))  max(Vertex(:,3))]);
          axis([min(min(Vertex))  max(max(Vertex))  min(min(Vertex))  max(max(Vertex))  min(min(Vertex))  max(max(Vertex))]);

        title(['tt:' num2str(tt,'%.2e'),' of ' num2str(num_step,'%.e'), '; dt:' num2str(dt,'%.e'),'; F_{rep}:' num2str(k),'; F_{adh}:' num2str(k_adh),'; \rho_{RSph}:' num2str(rho),'; speed:' num2str(v0),'; N:' num2str(num_part)]);  
%         xlabh = get(gca,'XLabel');
%         set(xlabh,'Position',get(xlabh,'Position'))% + [2 2 0])
%         ylabh = get(gca,'YLabel');
%         set(ylabh,'Position',get(ylabh,'Position'))% + [2.5 5 0])
%         xlabel('x') % x-axis label
%         ylabel('y') % y-axis label
%         zlabel('z')
        
        %quiver3(r(:,1),r(:,2),r(:,3),scale*r_dot(:,1),scale*r_dot(:,2),scale*r_dot(:,3),0,'MaxHeadSize', .8,'color','r');
        %quiver3(r(:,1),r(:,2),r(:,3),scale*F_track(:,1),scale*F_track(:,2),scale*F_track(:,3),0,'MaxHeadSize', .8,'color','r');
%         hold on

         
        %orientation of the particle (O); black:
        quiver3(r(:,1),r(:,2),r(:,3),scale*n(:,1),scale*n(:,2),scale*n(:,3),0,'MaxHeadSize', .8,'color','k');
        
%         %previous orentation of particles
%         quiver3(r_prev(:,1),r(:,2),r(:,3),scale*n_prev(:,1),scale*n_prev(:,2),scale*n_prev(:,3),0,'MaxHeadSize', .8,'color','g');
        %normalized velocity vector
        quiver3(r(:,1),r(:,2),r(:,3),scale*nr_dot(:,1),scale*nr_dot(:,2),scale*nr_dot(:,3),0,'MaxHeadSize', .8,'color','r');
        
%         %Direction of Force
%         quiver3(r(:,1),r(:,2),r(:,3),scale*F_track(:,1),scale*F_track(:,2),scale*F_track(:,3),0,'MaxHeadSize', .8,'color','y');

          %normal vector (N); blue:
%         quiver3(r(:,1),r(:,2),r(:,3),scale*Norm_vect(:,1),scale*Norm_vect(:,2),scale*Norm_vect(:,3),0,'MaxHeadSize', .8,'color','b');
        hold on
        %Binormal vector
        quiver3(r(:,1),r(:,2),r(:,3),scale*nr_dot_cross(:,1),scale*nr_dot_cross(:,2),scale*nr_dot_cross(:,3),0,'MaxHeadSize', .8,'color','b')
        hold on
     
        %Plot particle edges 
        plot3(r(:,1),r(:,2),r(:,3),'o','MarkerSize', 6,'MarkerEdgeColor','k');
        
        %-----Delauney Mesh to represent connectivity
%         hold on
%         tri2 = delaunayTriangulation(r(:,1),r(:,2),r(:,3));
%         [tri,Xb]= freeBoundary(tri2); 
%         h2=trimesh(tri,Xb(:,1),Xb(:,2),Xb(:,3),'FaceColor','none', 'EdgeColor', 'k','LineWidth', 2,'EdgeAlpha', 0.3);
%         hold off
        %pause
        
        axis vis3d;
%         view([-45 20]);
                view([0 90]);
%          camzoom(0.8);

        
        subplot(3,5,[9 15]); %EQUIRECTANGULAR Map PROJECTION - particle representation
        %-----------------------------------------------------------------------------
        
      
        %Definition of contourplot to show particle density
        F = TriScatteredInterp(thet1,phi1,number_neighbours);
        t_thet1= -pi:2*pi/100:pi;
%         t_phi1= -pi/2:pi/100:pi/2;
        t_phi1= -1.45:pi/100:1.45;
        [q_thet1,q_phi1]=meshgrid(t_thet1,t_phi1);
        q_neighbour=F(q_thet1,q_phi1);
        q_neighbour(q_neighbour==0)=nan;
        
        
        
        %plot contour according to particle density
        contourf(1.05*q_thet1,1.05*q_phi1, q_neighbour,20,'LineColor','none');hold on
        
        %adding color to white background
%         [whatever, conthandle] = contourf(1.05*q_thet1,1.05*q_phi1, q_neighbour,100,'LineColor','none');
%         
%         thismap = get(gcf, 'Colormap');
%         fillcol = thismap(1,:);
%         
%         set(gca, 'Color', fillcol);
%         
%         cc = get(conthandle, 'Children');
%         for K = 1:length(cc)
%             thiscol = get(cc(K), 'FaceColor');
%             if all(thiscol == 1) %all white and not string 'flat'
%                 set(cc(K), 'FaceColor', fillcol);
%             end
%         end
%         
%         
%         hold on
        
        
        %plot particle positions
        plot(thet1,phi1,'.','MarkerSize', 2,'Color','k');hold off
        
        hc2=colorbar;
        ylabel(hc2, 'number of nearest neighbours')
        %colorbar axis range
        caxis([3.5 8.5]);
        
        %Print title...
%         title('Particle density map','FontSize',10);
        axis([-pi  pi  -pi/2  pi/2]);

        %Axis label for equirectangular projection
        xlabel('\Theta (azimuth) -pi < x < pi') % x-axis label
        ylabel('\Phi (altitude) -pi/2 < z < pi/2') % y-axis label

        %prolate spheroid cylinder projection
%         xlabel('Theta (azimuth) -pi < x < pi') % x-axis label
%         ylabel('Z (elevation) min(Z) < Z < max(Z)') % y-axis label
        
        %oblate spheroid cylinder projection
%         ylabel('Phi (altitude) -pi/2 < y < pi/2') % x-axis label
%         xlabel('Z (elevation) min(Z) < Z < max(Z)') % y-axis label

%         hold off
        
        ax=subplot(3,5,[4 5]); %ORDER PARAMETER
        %------------------------------------------------------------
%                 en_v(tt/plotstep)=norm(sum(r_dot));

        
        %squared norm of each element in r_dot
        ner_dot=[];
        for i=1:num_part
        ner_dot=[ner_dot; norm(r_dot(i,:))^2];
        end
        
        %total energy in the system
         if 0.5*sum(ner_dot)<200
        en_v(tt/plotstep)=0.5.*sum(ner_dot);
         end
        
        
%         plot(v_order);

        x_en_v=0:length(en_v)-1;
        x_v_order=0:length(v_order)-1;

        plotyy(0:length(v_order)-1,v_order,0:length(en_v)-1,en_v)
        [hAx,hLine1,hLine2] = plotyy(x_v_order,v_order,x_en_v,en_v);
        

        title(['Mean OP: ', num2str(mean(v_order(1:tt/plotstep)),'%.2f'),'; \tau:', num2str(tau),'; R_{0}:', num2str(r_adh),'; R_{eq}:', num2str(2*sigma_n,'%.2f'), '\Phi_{Pack}:' num2str(phi_pack,'%.2f'), '  SurfArea: ' num2str(Area,'%.2f')],'FontSize',10);
        xlabel('timestep tt/10')

        ylabel(hAx(1),'order parameter') % left y-axis
        ylabel(hAx(2),'total energy') % right y-axis  
%           set(hAx(1),'ylim',[0 1])
%           set(hAx(2),'ylim',[0 10])   
%         set(ax(1),'Box','off')
%         set(ax(1),'YTick',[0:.5:1])
%         set(ax(2),'YTick',[0:.5:10])
        
%         grid on   
        
%         subplot(3,5,5); %Energy
        %------------------------------------------------------------
        
        
%         en_v(tt/plotstep)=norm(sum(r_dot));
        
%         tot_F(tt/plotstep)=norm(sum(F_track));
        
%         plot(en_v);
%         title(['\Phi_{Pack}:' num2str(phi_pack,'%.2f'), '  SurfArea: ' num2str(Area,'%.2f')],'FontSize',10);
%         xlabel('timestep tt/10')
%         ylabel('total energy/forces')
%         grid on 
%         
%         hold on
%         plot(tot_F, 'r');
        
%         x_en_v=0:length(en_v)-1;
%         x_tot_F=0:length(tot_F)-1;
% 
%         plotyy(0:length(en_v)-1,en_v,0:length(tot_F)-1,tot_F)
%         [hAx,hLine1,hLine2] = plotyy(x_en_v,en_v,x_tot_F,tot_F);
% 
%         title(['\Phi_{Pack}:' num2str(phi_pack,'%.2f'), '  SurfArea: ' num2str(Area,'%.2f')],'FontSize',10);
%         xlabel('timestep tt/10')
% 
%         ylabel(hAx(1),'total energy') % left y-axis
%         ylabel(hAx(2),'sum forces') % right y-axis
        
        
%         if tt>200
%         axis([0 num_step/plotstep 0 max(en_v(10:end))])
%         end
       
% r = gather(r);
% n=gather(n);
        file_name_plot = fullfile(folder_plots,['image',num2str(tt/plotstep),'.tiff']);
        if movie_making == 1

%             drawnow;       % updates figure
                    
            frame = getframe(hFig);
            writeVideo(writerObj,frame);
        else
            saveas(hFig,file_name_plot)
        end
%         close
%         if tt>1 && tt<num_step
% %         delete(hFig);
%         close
%         end
    end
end

if movie_making == 1
    close(writerObj);
    name_data(name_data == '.') = '_';
    saveas(hFig,fullfile(UpFolder,'movies',[name_data,'.fig']))
%     saveas(hFig,file_name_plot)

end
file_struct = fullfile(folder_particle_simula,[name_data,'.mat']);
%file_structCSV = fullfile(folder_particle_simula,[namestructure,'.csv']);
save(file_struct,'-struct','particle_info','-v7.3')
%struct2csv(file_struct,file_structCSV)
display(['Times elapsed in hour: ',num2str(toc/3600)])