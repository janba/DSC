#version 150

uniform mat4 MVMatrix;
uniform mat4 MVPMatrix;
uniform mat4 NormalMatrix;

uniform vec3 lightPos;

in vec3 position;
in vec3 normal;

out vec4 colourV;

void main()
{
    vec4 p = MVMatrix * vec4(position.xyz, 1.);
    vec3 N = normalize(mat3(NormalMatrix) * normal);
    vec3 L = normalize(lightPos - p.xyz);
    
    float s = max(dot(-N,L), 0.0);
    colourV = vec4(s, s, s, 1.);
    
    // Calculate position
    gl_Position = MVPMatrix * vec4(position.xyz, 1.);
}
