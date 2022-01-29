#include "SPH.h"

//structure to store positions and time
struct vec3t{
    public :
        vcl::vec3 p; // position
        float t;     // time
};

//keyframes is a class that stores keyframes and has methods to calculate position with time. It also contains a timer
class keyframes{
    public:

    //create default object
    keyframes();
    
    //initialize object with positions of keyframes and speed of objects
    void init(vcl::buffer<vec3> keyframe_position, float speed);

    //update position of object
    void update(float simulation_dtime);

    //retrieve and set arbitrary times for keyframes
    void configure_time(std::vector<float> times);
    std::vector<float> retrieve_time();

    //attributes
    vcl::buffer<vec3t> keyframesBuffer;
    vcl::timer_interval timer_keyframe;
    vcl::vec3 p;
    vcl::vec3 pp;
    float K_for_spline_cardinal;

    private:
    //private attributes
    int picked_object;
};

class keyframes_matrices{
    public:

    //create default object
    keyframes_matrices();

    void init(vcl::buffer<mat3> keyframes_matrice, float speed);
    void update(float simulation_dtime);

    mat3 matrix;

    private:
    keyframes column0;
    keyframes column1;
    keyframes column2;

};

//return index of keframes array for value t
size_t index_at_value(float t, vcl::buffer<vec3t> const& v);

//interpolate cardinal splines
vec3 cardinal_spline_interpolation(float t, float t0, float t1, float t2, float t3, const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, float K);
vec3 cardinal_spline_interpolation_deriv(float t, float t0, float t1, float t2, float t3, const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, float K);
