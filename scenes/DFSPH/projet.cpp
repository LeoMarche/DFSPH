
#include "projet.hpp"
#include "generate_terrain.hpp"

#ifdef SCENE_PROJET

using namespace vcl;

void scene_model::setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{
    glEnable(GL_DEPTH_TEST);

    std::vector<float> times = {-2.0f, 0.0f, 2.0f, 2.78f, 6.62f, 8.36f, 10.86f, 12.86f};
    vcl::buffer<vec3> keyframe_position_camera = {{-1.54f, -5.0f, -11.76f},{-1.54f, -5.0f, -11.76f},{-1.54f, -5.0f, -11.76f},{-4.4f,-3.0f,-10.17f},{2.19f,-0.89f,-7.71f},{3.48f,-0.94f,-4.74f},{19.89, 16.55f, 22.74f},{11.89, -5.0f, 1.74f},{-1.54f, -5.0f, -11.76f}};
    vcl::buffer<vec3> keyframe_rotation_camera= {{M_PI/2, -M_PI, 0.0f},{M_PI/2, -M_PI, 0.0f},{M_PI/2, -M_PI, 0.0f},{0.79f,-1.53f,-0.45f},{2.88f,-0.86f,3.05f-2*M_PI},{2.67f,-0.53f,-3.12f},{-2.57f+2*M_PI, -0.86f, -2.91},{M_PI/2, -0.86f, 0.0f},{M_PI/2, -M_PI, 0.0f}};

    camera_position_keyframes.init(keyframe_position_camera, 1.9f);
    camera_orientation_keyframes.init(keyframe_rotation_camera, 1.0f);
    camera_position_keyframes.configure_time(times);
    camera_orientation_keyframes.configure_time(camera_position_keyframes.retrieve_time());

    //initialize bird and keyframes
    vcl::buffer<vec3> keyframe_position = {{-12.0721f,29.6224f,-0.15637f}, {-10.1453f,26.7269f,10.1432f}, {-4.28616f,27.4374f,13.1684f}, {1.94839f,26.8544f,10.1764f}, {4.03007f,26.2481f,10.7397f}, {6.10103f,27.5563f,9.46138f}, {9.1524f,28.9992f,8.12696f}, {13.0638f,23.8273f,-2.26078f}, {8.45956f,25.732f,-8.97049f}, {0.573354f,27.6128f,-9.58968f}, {-6.28622f,28.4363f,-8.86022f}, {-10.8069f,28.6868f,-5.85482f}, {-12.0721f,29.6224f,-0.15637f}, {-10.1453f,26.7269f,10.1432f}, {-4.28616f,27.4374f,13.1684f}};
    bird.init(shaders, keyframe_position);   
    
    //initialize sph simulation
    sph.set_wind(wind);
    sph.set_gravity(gravity);
    sph.setStaticSphere(0.f, 0.f, 0.f, 0.f);
    sph.set_nr_of_particles(N_particles);
    sph.set_max_dens_iter(max_iter_div);
    sph.set_max_div_iter(max_iter_dens);
    sph.set_divergence_error(div_error);
    sph.set_density_error(dens_error);
    sph.set_timestep(time_factor);
    sph.reset();
    sph.init_positions(offset_pos, rows_cols[0], rows_cols[1]);

    // Display elements
    // ******************************************* //
    //Sphere for particles
    sphere = mesh_primitive_sphere(1.0f, {0,0,0}, 10UL, 10UL);
    sphere.uniform.color = {0,0,1.0f};
    sphere.uniform.transform.scaling = 0.02f;
    sphere.uniform.shading.diffuse = 1.5f;

    //skybox
    skybox = create_skybox();
    skybox.uniform.transform.rotation = rotation_from_axis_angle_mat3({1,0,0},-3.1415926f/2.0f);
    skybox.uniform.transform.translation = {-100,-100,100};
    skybox.uniform.transform.scaling =200;
    skybox.texture_id = create_texture_gpu(image_load_png("scenes/DFSPH/assets/skybox.png"));

    //terrain
    terrain = create_terrain();
    terrain.uniform.transform.translation = {0,0,0};
    terrain.uniform.shading.specular = 0.0f;
    terrain.uniform.shading.diffuse = 1.5f;
    terrain.texture_id = create_texture_gpu(image_load_png("scenes/DFSPH/assets/final_4k_rocks.png"));

    //trees
    float foliage_radius = 2.0f;
    tree_positions = {vec2(0.9f,0.9f)};
    add_random_pos(tree_positions, foliage_radius, 20, vec2(0.04f,0.58f), vec2(0.58f,0.96f));
    add_random_pos(tree_positions, foliage_radius, 8, vec2(0.8f, 0.95f), vec2(0.04f, 0.98f));

    trunc = create_cylinder(0.3,2.0);
    trunc.uniform.color = {0.35f,0.15f,0.0f};
    trunc.uniform.shading.specular = 0.05f;
    trunc.uniform.shading.diffuse = 1.5f;

    foliage = create_tree_foliage(foliage_radius,2.0f,1.0f);
    foliage.uniform.color = {0.1f,1.0f,0.1f};
    foliage.uniform.shading.specular = 0.05f;
    foliage.uniform.shading.diffuse = 1.5f;

    //shrooms
    shrooms_positions = {vec2(0.8f,0.8f)};
    foot_shrooms = create_cylinder(0.03f, foliage_radius/10);
    foot_shrooms.uniform.color = {1.0f, 1.0f, 1.0f};
    foot_shrooms.uniform.shading.diffuse=1.5f;
    foot_shrooms.uniform.shading.specular=0.05f;

    hat_shrooms = create_cone(0.2f, 0.2f, 0.1f);
    hat_shrooms.uniform.color={1.0f,0.0f,0.0f};
    hat_shrooms.uniform.shading.diffuse=1.5f;
    hat_shrooms.uniform.shading.specular=0.05f;

    add_random_pos(shrooms_positions, foliage_radius/10, 50, vec2(0.01f,0.58f), vec2(0.6f,0.99f));
    add_random_pos(shrooms_positions, foliage_radius/10, 10, vec2(0.8f, 0.98f), vec2(0.02f, 0.98f));

    //grass
    grass_positions = {vec2(0.02f,0.95f)};
    grass = mesh_load_file_obj("scenes/DFSPH/assets/Grass_Plane_001.obj");
    grass.texture_id = create_texture_gpu(image_load_png("scenes/DFSPH/assets/Plane_1_D.png"));
    grass.uniform.transform.scaling = 0.01f;
    grass.uniform.shading.specular = 0.0f;
    add_random_pos(grass_positions, 0.0f, 1000, vec2(0.02f,0.58f), vec2(0.6f,0.98f));
    add_random_pos(grass_positions, 0.0f, 80, vec2(0.8f, 0.98f), vec2(0.02f, 0.98f));

    //rocks
    rocks_positions = {vec2(0.85f, 0.85f)};
    rock = mesh_load_file_obj("scenes/DFSPH/assets/Rock_1.OBJ");
    rock.texture_id = create_texture_gpu(image_load_png("scenes/DFSPH/assets/Rock_1_Base_Color.png"));
    rock.uniform.transform.scaling = 0.05f;
    rock.uniform.shading.specular=0.0f;
    add_random_pos(rocks_positions, 0.5f, 20, vec2(0.02f,0.58f), vec2(0.6f,0.98f));
    add_random_pos(rocks_positions, 0.5f, 8, vec2(0.84f, 0.98f), vec2(0.02f, 0.98f));

    rocks2_positions = {vec2(0.5f,0.7f)};
    rock2 = mesh_load_file_obj("scenes/DFSPH/assets/rock1.obj");
    rock2.texture_id = create_texture_gpu(image_load_png("scenes/DFSPH/assets/rock1_diffuse2.png"));
    rock2.uniform.transform.scaling = 0.05f;
    rock2.uniform.shading.specular=0.0f;
    add_random_pos(rocks2_positions, 0.5f, 15, vec2(0.02f,0.58f), vec2(0.6f,0.98f));
    add_random_pos(rocks2_positions, 0.5f, 6, vec2(0.82f, 0.98f), vec2(0.02f, 0.98f));

    //initialize camera
    scene.camera.translation = keyframe_position_camera[0];
    scene.camera.orientation = rotation_from_axis_angle_mat3({0,0,1},keyframe_rotation_camera[0].z)*rotation_from_axis_angle_mat3({0,1,0},keyframe_rotation_camera[0].y)*rotation_from_axis_angle_mat3({1,0,0},keyframe_rotation_camera[0].x);
}





void scene_model::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{

    //get timestep
    float simulation_dtime = sph.get_timestep();

    //update camera
    camera_position_keyframes.update(simulation_dtime);
    camera_orientation_keyframes.update(simulation_dtime);
    scene.camera.translation = camera_position_keyframes.p;
    scene.camera.orientation = rotation_from_axis_angle_mat3({0,0,1}, camera_orientation_keyframes.p.z)*rotation_from_axis_angle_mat3({0,1,0}, camera_orientation_keyframes.p.y)*rotation_from_axis_angle_mat3({1,0,0}, camera_orientation_keyframes.p.x);
    
    //retrieve camera position (animation purpose)
    camera_position = scene.camera.translation;
    camera_rotation = scene.camera.orientation;

    set_gui();
    
    //draw skybox and terrain (comment/uncomment for textures)
    draw(skybox, scene.camera, shaders["mesh"]);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);

    //affichage d'une texture pour le terrain
    draw(terrain, scene.camera, shaders["mesh"]);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);//pas de texture pour les autres

    //trees
    for(unsigned int i=0; i<tree_positions.size(); i++){
        vec3 normale = -get_normale(tree_positions[i].x, tree_positions[i].y);
        mat3 rotation = rotation_between_vector_mat3({0,1,0}, normale);

        trunc.uniform.transform.translation = evaluate_terrain(tree_positions[i].x, tree_positions[i].y)-0.1f*normale;
        trunc.uniform.transform.rotation = rotation;
        foliage.uniform.transform.translation = trunc.uniform.transform.translation + 1.8f*normale;
        foliage.uniform.transform.rotation = rotation;

        draw(trunc, scene.camera, shaders["mesh"]);
        draw(foliage, scene.camera, shaders["mesh"]);
    }

    //shrooms
    for(unsigned int i=0; i<shrooms_positions.size(); i++){
        vec3 normale = -get_normale(shrooms_positions[i].x, shrooms_positions[i].y);
        mat3 rotation = rotation_between_vector_mat3({0,1,0}, normale);

        foot_shrooms.uniform.transform.translation = evaluate_terrain(shrooms_positions[i].x, shrooms_positions[i].y)-0.01f*normale;
        foot_shrooms.uniform.transform.rotation = rotation;
        hat_shrooms.uniform.transform.translation = foot_shrooms.uniform.transform.translation + 0.09f*normale;
        hat_shrooms.uniform.transform.rotation = rotation;

        draw(foot_shrooms, scene.camera, shaders["mesh"]);
        draw(hat_shrooms, scene.camera, shaders["mesh"]);
    }

    //rocks
    for(unsigned int i=0; i<rocks_positions.size(); i++){
        vec3 normale = -get_normale(rocks_positions[i].x, rocks_positions[i].y);
        mat3 rotation = rotation_between_vector_mat3({0,1,0}, normale);
        mat3 base_rotat = rotation_from_axis_angle_mat3({1,0,0},-M_PI/2);

        rock.uniform.transform.translation = evaluate_terrain(rocks_positions[i].x, rocks_positions[i].y)-0.001f*normale;
        rock.uniform.transform.rotation = rotation*rotation_from_axis_angle_mat3({0,1,0},float(i)/rocks_positions.size()*2*M_PI)*base_rotat;
        rock.uniform.transform.scaling = float(i)/rocks_positions.size()*0.05f;

        draw(rock, scene.camera, shaders["mesh"]);
    }
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);

    for(unsigned int i=0; i<rocks2_positions.size(); i++){
        vec3 normale = -get_normale(rocks2_positions[i].x, rocks2_positions[i].y);
        mat3 rotation = rotation_between_vector_mat3({0,1,0}, normale);
        mat3 base_rotat = rotation_from_axis_angle_mat3({1,0,0},-M_PI/2);

        rock2.uniform.transform.translation = evaluate_terrain(rocks2_positions[i].x, rocks2_positions[i].y)-0.001f*normale;
        rock2.uniform.transform.rotation = rotation*rotation_from_axis_angle_mat3({0,1,0},float(i)/rocks2_positions.size()*2*M_PI)*base_rotat;
        rock2.uniform.transform.scaling = float(i)/rocks2_positions.size()*0.05f;

        draw(rock2, scene.camera, shaders["mesh"]);
    }
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);

    //update and draw fluid particles
    sph.update();
    Float3* particle_pos = sph.get_particle_positions();
    for (int i = 0; i < sph.get_nr_of_particles() - 1; ++i) {
        vec3 particlePos = vec3(particle_pos->x[i],particle_pos->y[i],particle_pos->z[i]);
        sphere.uniform.transform.translation = particlePos;
        draw(sphere, scene.camera, shaders["mesh"]);
    }

    //update and draw birds
    bird.update(simulation_dtime);
    
    if(gui_scene.skeleton){
        bird.hierarchy_visual_debug.draw(bird.hierarchy, scene.camera);
    }
    
    if(gui_scene.surface){
        draw(bird.hierarchy, scene.camera);
    }
    if(gui_scene.surface){
        //draw(bird.target_position, scene.camera);
    }
    if(gui_scene.skeleton){
        bird.hierarchy_visual_debug.draw(bird.target_position, scene.camera);
    }
    if (is_first_frame){
        for (int i = 0;i<10;i++){
            bird.isFirstFrame(i);
            draw(bird.hierarchy, scene.camera);
        }
        is_first_frame = false;
    }
    else{
        for (int i=0;i<10;i++){
            bird.frame(i, simulation_dtime);
            draw(bird.hierarchy, scene.camera);
        }
        
    }

    //grass (png)
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(false);

    for(unsigned int i=0; i<grass_positions.size(); i++){
        vec3 normale = -get_normale(grass_positions[i].x, grass_positions[i].y);
        mat3 rotation = rotation_between_vector_mat3({0,1,0}, normale);

        grass.uniform.transform.translation = evaluate_terrain(grass_positions[i].x, grass_positions[i].y)-0.001f*normale;
        grass.uniform.transform.rotation = rotation*rotation_from_axis_angle_mat3({0,1,0},float(i)/grass_positions.size()*2*M_PI);
        grass.uniform.transform.scaling = float(i)/grass_positions.size()*0.005f;

        draw(grass, scene.camera, shaders["mesh"]);
    }
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    glDepthMask(true);

    //take screencapture of openGL and save it (uncomment to start)
    /*if(camera_position_keyframes.timer_keyframe.t-last_time>0.02f){
        save_screenshot_to_file("rendus/"+std::to_string(number_of_frame)+".tga",1920,1080);
        last_time = camera_position_keyframes.timer_keyframe.t;
        number_of_frame += 1;
    }*/
    if(number_of_frame==751){
        exit(0);
    }
}

void scene_model::set_gui(){

    //options for drawing
    ImGui::Checkbox("Surface", &gui_scene.surface);
    ImGui::Checkbox("Skeleton", &gui_scene.skeleton);

    //options for birds
    ImGui::SliderFloat("Time keyframe", &camera_position_keyframes.timer_keyframe.t, camera_position_keyframes.timer_keyframe.t_min, camera_position_keyframes.timer_keyframe.t_max);
    ImGui::SliderFloat("K", &bird.birdKeyframes.K_for_spline_cardinal, 0.0f, 1.0f);

    ImGui::Text("Display: "); ImGui::SameLine();
    ImGui::Checkbox("keyframe", &gui_scene.display_keyframe); ImGui::SameLine();
    ImGui::Checkbox("polygon", &gui_scene.display_polygon);

    //options to print keyframes or camera info in stdout
    if( ImGui::Button("Print Keyframe") )
    {
        std::cout<<"keyframe_position={";
        for(size_t k=0; k<bird.birdKeyframes.keyframesBuffer.size(); ++k)
        {
            const vec3& p = bird.birdKeyframes.keyframesBuffer[k].p;
            std::cout<< "{"<<p.x<<"f,"<<p.y<<"f,"<<p.z<<"f}";
            if(k<bird.birdKeyframes.keyframesBuffer.size()-1)
                std::cout<<", ";
        }
        std::cout<<"}"<<std::endl;
    }
    if( ImGui::Button("Print Camera Position") )
    {
        std::cout<<camera_position<<std::endl;
    }
    if( ImGui::Button("Print Camera Orientation") )
    {
        std::cout<<rotationMatrixToEulerAngles(camera_rotation)<<std::endl;
    }
    if( ImGui::Button("Print Time") )
    {
        std::cout<<camera_position_keyframes.timer_keyframe.t<<std::endl;
    }
}

//function that create skybox
mesh create_skybox()
{
    // Number of samples of the terrain is N x N 
    mesh skybox; // temporary skybox storage (CPU only)
    skybox.position.resize(14);
    skybox.texture_uv.resize(14);
    float e=0.001f;
    skybox.position  [0] = {0,0,1};
    skybox.texture_uv[0] = {0,1.0f/3.0f};
    skybox.position  [1] = {0,0,0};
    skybox.texture_uv[1] = {0,2.0f/3.0f};
    skybox.position  [2] = {0,0,1};
    skybox.texture_uv[2] = {1.0f/4.0f,0.0f/3.0f};
    skybox.position  [3] = {1,0,1};
    skybox.texture_uv[3] = {1.0f/4.0f,1.0f/3.0f};
    skybox.position  [4] = {1,0,0};
    skybox.texture_uv[4] = {1.0f/4.0f,2.0f/3.0f};
    skybox.position  [5] = {0,0,0};
    skybox.texture_uv[5] = {1.0f/4.0f,1-e};
    skybox.position  [6] = {0,1,1};
    skybox.texture_uv[6] = {2.0f/4.0f,0.0f/3.0f};
    skybox.position  [7] = {1,1,1};
    skybox.texture_uv[7] = {2.0f/4.0f,1.0f/3.0f};
    skybox.position  [8] = {1,1,0};
    skybox.texture_uv[8] = {2.0f/4.0f,2.0f/3.0f};
    skybox.position  [9] = {0,1,0};
    skybox.texture_uv[9] = {2.0f/4.0f,1-e};
    skybox.position  [10] = {0,1,1};
    skybox.texture_uv[10] = {3.0f/4.0f,1.0f/3.0f};
    skybox.position  [11] = {0,1,0};
    skybox.texture_uv[11] = {3.0f/4.0f,2.0f/3.0f};
    skybox.position  [12] = {0,0,1};
    skybox.texture_uv[12] = {1-e,1.0f/3.0f};
    skybox.position  [13] = {0,0,0};
    skybox.texture_uv[13] = {1-e,2.0f/3.0f};

    skybox.connectivity.push_back({1,0,3});
    skybox.connectivity.push_back({4,1,3});
    skybox.connectivity.push_back({3,2,6});
    skybox.connectivity.push_back({7,3,6});
    skybox.connectivity.push_back({4,3,7});
    skybox.connectivity.push_back({8,4,7});
    skybox.connectivity.push_back({5,4,8});
    skybox.connectivity.push_back({9,5,8});
    skybox.connectivity.push_back({8,7,10});
    skybox.connectivity.push_back({11,8,10});
    skybox.connectivity.push_back({11,10,12});
    skybox.connectivity.push_back({13,11,12});


    return skybox;
}

//function that generate terrain mesh
mesh create_terrain()
{
    // Number of samples of the terrain is N x N
    const size_t N = 1000;
 
    mesh terrain; // temporary terrain storage (CPU only)
    terrain.position.resize(N*N);
    terrain.texture_uv.resize(N*N);

    // Fill terrain geometry
    for(size_t ku=0; ku<N; ++ku)
    {
        for(size_t kv=0; kv<N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku/(N-1.0f);
            const float v = kv/(N-1.0f);

            // Compute coordinates
            terrain.position[kv+N*ku] = evaluate_terrain(u,v);
            terrain.texture_uv[kv+N*ku]={(float) ku/(N-1),(float) kv/(N-1)};
        }
    }
    // Generate triangle organization
    //  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
    const unsigned int Ns = N;
    for(unsigned int ku=0; ku<Ns-1; ++ku)
    {
        for(unsigned int kv=0; kv<Ns-1; ++kv)
        {
            const unsigned int idx = kv + N*ku; // current vertex offset

            const uint3 triangle_1 = {idx, idx+1+Ns, idx+1};
            const uint3 triangle_2 = {idx, idx+Ns, idx+1+Ns};

            terrain.connectivity.push_back(triangle_1);
            terrain.connectivity.push_back(triangle_2);

            
        }
    }
    return terrain;
}

mesh create_tree_foliage(float radius, float height, float z_offset)
{
    mesh m = create_cone(radius, height, 0);
    m.push_back( create_cone(0.7*radius, height, z_offset) );
    m.push_back( create_cone(0.5*radius, height, 2*z_offset) );

    return m;
}

mesh create_cone(float radius, float height, float z_offset)
{
    mesh m;

    // conical structure
    // *************************** //

    const size_t N = 20;

    // geometry
    for(size_t k=0; k<N; ++k)
    {
        const float u = k/float(N);
        const vec3 p = {radius*std::cos(2*3.14f*u), 0.0f, radius*std::sin(2*3.14f*u)};
        m.position.push_back( p+vec3{0,z_offset,0} );
    }
    // apex
    m.position.push_back({0,height+z_offset,0});

    // connectivity
    const unsigned int Ns = N;
    for(unsigned int k=0; k<Ns; ++k) {
        m.connectivity.push_back( {k , (k+1)%Ns, Ns} );
    }

    // close the bottom of the cone
    // *************************** //

    // Geometry
    for(size_t k=0; k<N; ++k)
    {
        const float u = k/float(N);
        const vec3 p = {radius*std::cos(2*3.14f*u), 0.0f, radius*std::sin(2*3.14f*u)};
        m.position.push_back( p+vec3{0,z_offset,0} );
    }
    // central position
    m.position.push_back( {0,z_offset,0} );

    // connectivity
    for(unsigned int k=0; k<Ns; ++k)
        m.connectivity.push_back( {k+Ns+1, (k+1)%Ns+Ns+1, 2*Ns+1} );

    return m;
}

mesh create_cylinder(float radius, float height)
{
    mesh m;

    // Number of samples
    const size_t N = 20;

    // Geometry
    for(size_t k=0; k<N; ++k)
    {
        const float u = k/float(N);
        const vec3 p = {radius*std::cos(2*3.14f*u), 0.0f, radius*std::sin(2*3.14f*u)};
        m.position.push_back( p );
        m.position.push_back( p+vec3(0,height,0) );
    }

    // Connectivity
    for(size_t k=0; k<N; ++k)
    {
        const unsigned int u00 = 2*k;
        const unsigned int u01 = (2*k+1)%(2*N);
        const unsigned int u10 = (2*(k+1))%(2*N);
        const unsigned int u11 = (2*(k+1)+1) % (2*N);

        const uint3 t1 = {u00, u10, u11};
        const uint3 t2 = {u00, u11, u01};
        m.connectivity.push_back(t1);
        m.connectivity.push_back(t2);
    }

    return m;
}

void add_random_pos(std::vector<vcl::vec2> &tab, float radius, int number, vec2 limits_x, vec2 limits_y){

    for(int i=0; i<number; i++){
        bool good = false;
        float u;
        float v;
        while(!good){
            u = rand_interval(limits_x.x, limits_x.y);
            v = rand_interval(limits_y.x, limits_y.y);
            for(unsigned int i=0; i<tab.size(); i++){
                good = true;
                if(50*sqrt((u-tab[i].x)*(u-tab[i].x)+(v-tab[i].y)*(v-tab[i].y))<2*radius){
                    good = false;
                    break;
                }
            }
        }
        tab.push_back(vec2(u,v));
    }
}

void save_screenshot_to_file(std::string filename, int windowWidth, int windowHeight) {    
    const int numberOfPixels = windowWidth * windowHeight * 3;
    unsigned char pixels[numberOfPixels];

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, windowWidth, windowHeight, GL_BGR, GL_UNSIGNED_BYTE, pixels);

    FILE *outputFile = fopen(filename.c_str(), "w");
    short header[] = {0, 2, 0, 0, 0, 0, (short) windowWidth, (short) windowHeight, 24};

    fwrite(&header, sizeof(header), 1, outputFile);
    fwrite(pixels, numberOfPixels, 1, outputFile);
    fclose(outputFile);

    printf("Finish writing to file %s \n");
    std::cout<<filename<<std::endl;
}


// Calculates rotation matrix to euler angles
vec3 rotationMatrixToEulerAngles(mat3 R)
{
    
    float sy = sqrt(R.xx * R.xx +  R.yx * R.yx );

    bool singular = sy < 1e-6; // If

    float x, y, z;
    if (!singular)
    {
        x = atan2(R.zy , R.zz);
        y = atan2(-R.zx, sy);
        z = atan2(R.yx, R.xx);
    }
    else
    {
        x = atan2(-R.yz, R.yy);
        y = atan2(-R.zx, sy);
        z = 0;
    }
    return vec3(x, y, z);  
}
#endif
