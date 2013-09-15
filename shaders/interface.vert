#version 150

uniform mat4 MVP;

in vec3 position;

out vec4 colourV;

void main (void)
{
    colourV = vec4(1.,0.,0.,1.);
    gl_Position = MVP*vec4(position.xyz, 1.);
}
