#include "keyframes.hpp"

keyframes::keyframes(){}

void keyframes::init(vcl::buffer<vec3> keyframe_position, float speed){

    //push all keyframes positions in object
    for (unsigned int i = 0;i<keyframe_position.size();i++){
            keyframesBuffer.push_back({keyframe_position[i],0});
        }

        float tsum=0.0f;
        keyframesBuffer[0].t=tsum;
        for (long unsigned int i = 1; i<keyframesBuffer.size();i++){
            float sum = (keyframesBuffer[i].p[0]-keyframesBuffer[i-1].p[0])*(keyframesBuffer[i].p[0]-keyframesBuffer[i-1].p[0]);
            sum += (keyframesBuffer[i].p[1]-keyframesBuffer[i-1].p[1])*(keyframesBuffer[i].p[1]-keyframesBuffer[i-1].p[1]);
            sum += (keyframesBuffer[i].p[2]-keyframesBuffer[i-1].p[2])*(keyframesBuffer[i].p[2]-keyframesBuffer[i-1].p[2]);
            tsum+= std::sqrt(sum)/speed;
            keyframesBuffer[i].t = tsum;
        }
    
    //initialize timer and attributes
    keyframesBuffer[0].t=keyframesBuffer[keyframesBuffer.size()-2].t-keyframesBuffer[keyframesBuffer.size()-1].t;
    timer_keyframe.t_min = keyframesBuffer[1].t;                   // first time of the keyframe
    timer_keyframe.t_max = keyframesBuffer[keyframesBuffer.size()-2].t;  // last time of the keyframe
    timer_keyframe.t = timer_keyframe.t_min;
    timer_keyframe.scale = 0;
    picked_object=-1;
    K_for_spline_cardinal =0.5f;//spline
}

void keyframes::update(float simulation_dtime){

    //update timer
    timer_keyframe.t+=simulation_dtime;
    timer_keyframe.update();

    const int idx = index_at_value(timer_keyframe.t, keyframesBuffer);

    // Preparation of data for the linear interpolation
    // Parameters used to compute the linear interpolation
    const float t0 = keyframesBuffer[idx-1].t; // t_{i-1}
    const float t1 = keyframesBuffer[idx  ].t; // = t_i
    const float t2 = keyframesBuffer[idx+1].t; // = t_{i+1}
    const float t3 = keyframesBuffer[idx+2].t; // = t_{i+2}

    const vec3& p0 = keyframesBuffer[idx-1].p; // = p_{i-}
    const vec3& p1 = keyframesBuffer[idx  ].p; // = p_i
    const vec3& p2 = keyframesBuffer[idx+1].p; // = p_{i+1}
    const vec3& p3 = keyframesBuffer[idx+2].p; // = p_{i+2}

    // Cardinal_spline_interpolation
    p = cardinal_spline_interpolation(timer_keyframe.t,t0,t1,t2,t3,p0,p1,p2,p3,K_for_spline_cardinal);
    pp = cardinal_spline_interpolation_deriv(timer_keyframe.t,t0,t1,t2,t3,p0,p1,p2,p3,K_for_spline_cardinal);
}

void keyframes::configure_time(std::vector<float> times){
    for(unsigned int i=0; i<times.size(); i++){
        if(i<keyframesBuffer.size()){
            keyframesBuffer[i].t = times[i];
        }
    }
    timer_keyframe.t_min = keyframesBuffer[1].t;
    timer_keyframe.t_max = keyframesBuffer[keyframesBuffer.size()-2].t;

}

std::vector<float> keyframes::retrieve_time(){
    std::vector<float> ret;
    for(unsigned int i=0; i<keyframesBuffer.size();i++){
        ret.push_back(keyframesBuffer[i].t);
    }
    return ret;
}

keyframes_matrices::keyframes_matrices(){}

void keyframes_matrices::init(vcl::buffer<mat3> buffer_keyframes, float speed){
    vcl::buffer<vec3> init_c0;
    vcl::buffer<vec3> init_c1;
    vcl::buffer<vec3> init_c2;

    for(unsigned int i = 0; i<buffer_keyframes.size(); i++){
        init_c0.push_back(buffer_keyframes[i].col(0));
        init_c1.push_back(buffer_keyframes[i].col(1));
        init_c2.push_back(buffer_keyframes[i].col(2));
    }

    //init columns
    column0.init(init_c0, speed);
    column1.init(init_c1, speed);
    column2.init(init_c2, speed);
}

void keyframes_matrices::update(float simulation_dtime){

    //update columns
    column0.update(simulation_dtime);
    column1.update(simulation_dtime);
    column2.update(simulation_dtime);

    //update matrix
    matrix = mat3(column0.pp, column1.pp, column2.pp);
}


size_t index_at_value(float t, vcl::buffer<vec3t> const& v)
{
    const size_t N = v.size();
    assert(N>=2);
    assert(t>=v[0].t);
    assert(t<v[N-1].t);
    size_t k=1;
    while( v[k+1].t<t ){
        ++k;
    }
    return k;
}

vec3 cardinal_spline_interpolation(float t, float t0, float t1, float t2, float t3, const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, float K){

    const float s = (t-t1)/(t2-t1);
    const vec3 d1 = 2*K*(p2-p0)/(t2-t0);
    const vec3 d2 = 2*K*(p3-p1)/(t3-t1);
    const vec3 p = (2*s*s*s-3*s*s+1)*p1+(s*s*s-2*s*s+s)*d1+(-2*s*s*s+3*s*s)*p2+(s*s*s-s*s)*d2;
    return p;
}

vec3 cardinal_spline_interpolation_deriv(float t, float t0, float t1, float t2, float t3, const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3, float K){

    const float s = (t-t1)/(t2-t1);
    const vec3 d1 = 2*K*(p2-p0)/(t2-t0);
    const vec3 d2 = 2*K*(p3-p1)/(t3-t1);
    const vec3 p = (1/(t2-t1))*(2*3*s*s-3*2*s)*p1+(1/(t2-t1))*(3*s*s-2*2*s+1)*d1+(1/(t2-t1))*(-2*3*s*s+3*2*s)*p2+(1/(t2-t1))*(3*s*s-2*s)*d2;
    return p;
}

