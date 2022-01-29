#include "keyframes.hpp"

class Bird{
    public:

    //create empry object
    Bird();

    //initialize object with shader and keyframes position for birdKeyframes
    void init(std::map<std::string, GLuint> shaders, vcl::buffer<vcl::vec3> keyframe_position);

    //function that updates individual positions (run only for first frame)
    void isFirstFrame(int i);

    //function that updates individual positions
    void frame(int i, float simulation_dtime);

    //update birds
    void update(float simulation_dtime);

    //attributes
    const float speed = 5.0f;
    vcl::hierarchy_mesh_drawable hierarchy;
    vcl::hierarchy_mesh_drawable target_position;
    vcl::buffer<vec3> actual_position;
    vcl::hierarchy_mesh_drawable_display_skeleton hierarchy_visual_debug;

    //keyframes of bird
    keyframes birdKeyframes;

    private:

    //rotation matrixes
    mat3 R_head;
    mat3 R_high_wing;
    mat3 R_low_wing;
    mat3 R_high_wing_r;
    mat3 R_low_wing_r;

    //size attibutes
    const float size_head = 0.8f;
    const float top_wing_length = 1.2f;
    const float middle_wing_length = 0.8f;
    const float top_wing_width = 1.6f;
    const float middle_wing_width = 1.3f;
    const float bottom_wing_width = 0.8f;

    //internal meshs
    mesh_drawable body2;
    mesh_drawable head2;
    mesh_drawable bec2;
    mesh_drawable eye2;
    mesh_drawable wing_top2;
    mesh_drawable wing_middle2;
    

};

//function to get rotation matrix from pp
mat3 pp_to_rotation_matrix(vec3 pp);

//function to get rotation matrix from pp only along y
mat3 pp_to_rotation_matrix_y_only(vec3 pp);