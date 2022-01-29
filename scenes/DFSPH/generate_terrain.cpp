#include "generate_terrain.hpp"

float evaluate_terrain_y(float u, float v)
{
    float func0 = 10.0f;
    float func1 = std::min(40.0f,250*std::abs(u-0.7f)+50*std::abs(v-0.5f));
    float func2 = std::min(250*((u-0.5)*(u-0.5)*(u-0.5)+0.5*(v-0.5)*(v-0.5f)+0.5*(u-0.5))-15,30.0);

    float height = 6.0f;
    float scaling = 1.5f;
    int octave = 7;
    float persistency = 0.4f;
    vec2 offset = vec2(0.1f, -0.25f);
    float noise = perlin((u-offset.x)*scaling,(v-offset.y)*scaling,octave,persistency,2.0F);
    //float noise = 0.01f;
    
    if(std::abs(u-0.5f)>0.5f || std::abs(v-0.5f)>0.5f){
        return func0 + (noise-1)*height;
    }
    else if(v-0.55>0){
        return (10*(v-0.5)*func1)/12.0 + (noise-1)*height;
    }
    else if (v-0.5>0){
        return (10*(v-0.5)*func1+2*std::abs(v-1)*func2*std::abs(v-0.55))/12.0 + (noise-1)*height;
    }
    return 2*std::abs(v-1)*func2*std::abs(v-0.55)/12.0+std::abs(v-0.5f)*func2 + (noise-1)*height;

}

vec3 get_normale(float u, float v){
    float grad_discrete_step = 0.001f;
	float y1 = evaluate_terrain_y(u-grad_discrete_step,v);
	float y2 = evaluate_terrain_y(u+grad_discrete_step,v);
	float y3 = evaluate_terrain_y(u,v-grad_discrete_step);
	float y4 = evaluate_terrain_y(u,v+grad_discrete_step);
	vec3 normale = cross(vec3(grad_discrete_step*100,y4-y3,0),vec3(0,y2-y1,grad_discrete_step*100));
    normale = normalize(normale);
    return normale;
}

// Evaluate 3D position of the terrain for any (u,v) \in [0,1]
vec3 evaluate_terrain(float u, float v)
{
    const float z = 50*(u-0.5f);
    const float x = 50*(v-0.5f);
    float y = evaluate_terrain_y(u,v);
    return {x,y,z};
}

