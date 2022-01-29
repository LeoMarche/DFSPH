#include "bird.hpp"

Bird::Bird(){}

void Bird::init(std::map<std::string, GLuint> shaders, vcl::buffer<vcl::vec3> keyframe_position){

    //initialize meshs
    mesh_drawable body2 = mesh_drawable( mesh_primitive_sphere(1, {0,0,0}, 40, 40));
    body2.uniform.transform.scaling_axis = {1.2f,1,2};

    mesh_drawable head2 = mesh_drawable(mesh_primitive_sphere(size_head, {0,0,0}, 40, 40));

    mesh_drawable bec2 = mesh_drawable(mesh_primitive_cone(size_head*0.4,{0,0,0},size_head*vec3(0,0,0.8f)));
    bec2.uniform.color = {1.0f,0.5f,0};//orange

    mesh_drawable eye2 = mesh_drawable(mesh_primitive_sphere(size_head*0.1f, {0,0,0}, 40, 40));
    eye2.uniform.color = {0,0,0};//black

    mesh_drawable wing_top2 = mesh_drawable(mesh_primitive_quad({0,0,-0.5f*top_wing_width},{0,0,0.5f*top_wing_width},{top_wing_length,0,0.5f*middle_wing_width},{top_wing_length ,0,-0.5f*middle_wing_width}));

    mesh_drawable wing_middle2 = mesh_drawable(mesh_primitive_quad({0,0,-0.5f*middle_wing_width},{0,0,0.5f*middle_wing_width},{middle_wing_length,0,0.5f*bottom_wing_width},{middle_wing_length ,0,-0.5f*bottom_wing_width}));
    
    //initialize hierarchy
    hierarchy.add(body2,"body");
    hierarchy.add(head2,"head","body",vec3(0.0f,0.8f,2.0f));
    hierarchy.add(bec2,"bec","head",size_head*vec3(0,0,0.9f));
    hierarchy.add(eye2, "eye_left", "head" , size_head * 1.1f * vec3( 1/3.0f, 1/2.0f, 1/1.5f));
    hierarchy.add(eye2, "eye_right", "head", size_head * 1.1f * vec3(-1/3.0f, 1/2.0f, 1/1.5f));
    hierarchy.add(wing_top2,"wing_high_left","body",1.2f*0.8f*vec3(1,0,0));
    hierarchy.add(wing_middle2,"wing_low_left","wing_high_left",top_wing_length*vec3(1,0,0));
    
    wing_top2.uniform.transform.scaling_axis = {-1,1,1};
    wing_middle2.uniform.transform.scaling_axis = {-1,1,1};

    hierarchy.add(wing_top2,"wing_high_right","body",-1.2f*0.8f*vec3(1,0,0));
    hierarchy.add(wing_middle2,"wing_low_right","wing_high_right",-top_wing_length*vec3(1,0,0));
    hierarchy["body"].transform.scaling = 0.4f;
    hierarchy.set_shader_for_all_elements(shaders["mesh"]);

   //initialize target_position
    target_position.add(head2,"center");
    for (int i =0;i<5;i++){
        target_position.add(head2,to_string(2*i),"center",vec3(1.0f*(i+1),0.0f,-2.0f*(i+1)));
        target_position.add(head2,to_string(2*i+1),"center",vec3(-1.0f*(i+1),0.0f,-2.0f*(i+1)));
    }
    target_position.set_shader_for_all_elements(shaders["mesh"]);
    hierarchy_visual_debug.init(shaders["segment_im"], shaders["mesh"]);

    //initialize keyframes for bird
    birdKeyframes.init(keyframe_position, speed);
}

void Bird::update(float simulation_dtime){

    //update keyframes
    birdKeyframes.update(simulation_dtime);
    float t = birdKeyframes.timer_keyframe.t;

    //update rotation and hierarchy for bird
    R_head = rotation_from_axis_angle_mat3({1,0,0}, 0.2*std::sin(2*3.14*t));
    R_high_wing = rotation_from_axis_angle_mat3({0,0,1},1.2f*std::sin(4*3.14f*(t)));
    R_low_wing = rotation_from_axis_angle_mat3({0,0,1},0.6f*std::sin(4*3.14f*(t)));
    R_high_wing_r = rotation_from_axis_angle_mat3({0,0,1},1.2f*std::sin(4*3.14f*(t-0.25f)));
    R_low_wing_r = rotation_from_axis_angle_mat3({0,0,1},0.6f*std::sin(4*3.14f*(t-0.25f)));
    hierarchy["body"].transform.translation = birdKeyframes.p;
    hierarchy["body"].transform.rotation = pp_to_rotation_matrix(birdKeyframes.pp);
    hierarchy["head"].transform.rotation = R_head;
    hierarchy["wing_high_left"].transform.rotation = R_high_wing;
    hierarchy["wing_low_left"].transform.rotation = R_low_wing;
    hierarchy["wing_high_right"].transform.rotation = R_high_wing_r;
    hierarchy["wing_low_right"].transform.rotation = R_low_wing_r;
    hierarchy.update_local_to_global_coordinates();

    //update target_position
    target_position["center"].transform.translation = birdKeyframes.p;
    target_position["center"].transform.rotation = pp_to_rotation_matrix_y_only(birdKeyframes.pp);
    target_position.update_local_to_global_coordinates();
}

void Bird::isFirstFrame(int i){
    actual_position.push_back(target_position[to_string(i)].global_transform.translation);
    hierarchy["body"].transform.translation = actual_position[i];
    hierarchy.update_local_to_global_coordinates();
}

void Bird::frame(int i, float simulation_dtime){
    vec3 target = target_position[to_string(i)].global_transform.translation;
    vec3 act = actual_position[i];
    float distance = (act.x-target.x)*(act.x-target.x)+(act.y-target.y)*(act.y-target.y)+(act.z-target.z)*(act.z-target.z);
    distance = std::sqrt(distance);
        if (speed*simulation_dtime>distance){
            actual_position[i]=target;
        }
        else{
            actual_position[i]+=(std::min(1.7f*speed*simulation_dtime,distance))*(target-act)/distance;
        }
        hierarchy["body"].transform.translation = actual_position[i];
        hierarchy["body"].transform.rotation = pp_to_rotation_matrix(target-act);
        hierarchy.update_local_to_global_coordinates();
}

mat3 pp_to_rotation_matrix(vec3 pp){
    if (pp.z==0&&pp.x==0)return rotation_from_axis_angle_mat3({0,0,1},-3.14f/2.0f);
    int signe = (pp.x>0) ? 1 : -1;
    mat3 rotationy = rotation_from_axis_angle_mat3({0,1,0}, signe *std::acos(pp.z/std::sqrt(pp.x*pp.x+pp.z*pp.z)));
    mat3 rotationx = rotation_from_axis_angle_mat3({1,0,0}, -std::atan(pp.y/std::sqrt(pp.x*pp.x+pp.z*pp.z)));

    return rotationy*rotationx;
}

mat3 pp_to_rotation_matrix_y_only(vec3 pp){
    if (pp.z==0&&pp.x==0)return rotation_from_axis_angle_mat3({0,0,1},0);
    int signe = (pp.x>0) ? 1 : -1;
    mat3 rotationy = rotation_from_axis_angle_mat3({0,1,0}, signe *std::acos(pp.z/std::sqrt(pp.x*pp.x+pp.z*pp.z)));

    return rotationy;
}