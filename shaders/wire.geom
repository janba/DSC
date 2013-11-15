#version 150

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

in vec4 colourV[3];

out vec4 colour;
out vec3 dist;

const int wscale = 512;

void main()
{
    vec4 p0 = wscale * gl_in[0].gl_Position/gl_in[0].gl_Position.w;
    vec4 p1 = wscale * gl_in[1].gl_Position/gl_in[1].gl_Position.w;
    vec4 p2 = wscale * gl_in[2].gl_Position/gl_in[2].gl_Position.w;
    
	vec2 v0 = p2.xy - p1.xy;
	vec2 v1 = p2.xy - p0.xy;
	vec2 v2 = p1.xy - p0.xy;
	float area = abs(v1.x*v2.y - v1.y * v2.x);
    
	dist = vec3(area/length(v0),0,0);
	gl_Position = gl_in[0].gl_Position;
	colour = colourV[0];
	EmitVertex();
    
	dist = vec3(0,area/length(v1),0);
	gl_Position = gl_in[1].gl_Position;
	colour = colourV[1];
	EmitVertex();
    
	dist = vec3(0,0,area/length(v2));
	gl_Position = gl_in[2].gl_Position;
	colour = colourV[2];
	EmitVertex();
    
	EndPrimitive();
}
