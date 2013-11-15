#version 150

uniform vec4 wireCol;

in vec3 dist;
in vec4 colour;

out vec4 fragColour;

void main(void)
{
    float d = min(dist.x,min(dist.y,dist.z));
	float I = exp2(-0.5*d*d);
	fragColour = I*wireCol + (1.0 - I)*colour;
}
