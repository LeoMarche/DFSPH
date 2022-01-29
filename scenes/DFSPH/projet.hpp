#pragma once

#include "SPH.h"
#include "bird.hpp"

#ifdef SCENE_PROJET

//function to create skybox
mesh create_skybox();

//functions to create terrain
mesh create_terrain();
vec3 evaluate_terrain(float u, float v);
float evaluate_terrain_y(float u, float v);
mesh create_tree_foliage(float radius, float height, float z_offset);
mesh create_cone(float radius, float height, float z_offset);
mesh create_cylinder(float radius, float height);
void add_random_pos(std::vector<vcl::vec2> &tab, float radius, int number, vec2 limits_x, vec2 limits_y);

//function to save screenshot of openGL window to file
void save_screenshot_to_file(std::string filename, int windowWidth, int windowHeight);

//function for retrieve angles
float norm(mat3 m);
bool isRotationMatrix(mat3 R);
vec3 rotationMatrixToEulerAngles(mat3 R);

//structure for gui
struct gui_scene_structure
{
    bool wireframe;
    bool surface     = true;
    bool skeleton    = false;
    bool display_keyframe = false;
    bool display_polygon  = false;
};

//structure for scene
struct scene_model : scene_base
{

    //void to setup elements and draw frames
    void setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);

    //number of particles
    int N_particles=1000;

    //Object for particles simulation
    SPH sph;
    Float3s wind{ 0.f,0.f,0.f };
    const float gravity=9.81f;
    const int max_iter_div = 200, max_iter_dens = 200;
    const float dens_error = 0.1f, div_error = 0.01f;
	const float max_error = 0.2f, min_error = 0.01f;
	const float time_factor = 0.4f;
    Float3s offset_pos = { 20.0f,10.0f,9.5f };
	const int rows_cols[2] = { 40, 156 };
    vcl::mesh_drawable sphere;

    //Object for birds flying
    Bird bird;
    bool is_first_frame =true;

    //keyframes for camera
    keyframes camera_position_keyframes;
    keyframes camera_orientation_keyframes;

    //Terrain & skybox
    vcl::mesh_drawable skybox;

    vcl::mesh_drawable terrain;

    vcl::mesh_drawable trunc;
    vcl::mesh_drawable foliage;
    vcl::mesh_drawable foot_shrooms;
    vcl::mesh_drawable hat_shrooms;
    vcl::mesh_drawable grass;

    vcl::mesh_drawable rock;
    vcl::mesh_drawable rock2;


    std::vector<vcl::vec2> tree_positions;
    std::vector<vcl::vec2> shrooms_positions;
    std::vector<vcl::vec2> rocks_positions;
    std::vector<vcl::vec2> rocks2_positions;
    std::vector<vcl::vec2> grass_positions;


    //vec with position of camera (for animation)
    vcl::vec3 camera_position;
    vcl::mat3 camera_rotation;

    //gui elements
    void set_gui();
    gui_scene_structure gui_scene;

    //simulation time used to store pictures
    float number_of_frame = 0;  
    float last_time = 0.0f;

};






#endif
