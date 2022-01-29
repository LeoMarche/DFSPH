#include "main/scene_base/base.hpp"
using namespace vcl;

float evaluate_terrain_y(float u, float v);
vec3 evaluate_terrain(float u, float v);
vec3 get_normale(float u, float v);
float Get2DPerlinNoiseValue(float x, float y, float res);