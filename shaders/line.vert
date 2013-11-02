#version 150

uniform mat4 MVMatrix;
uniform mat4 NormalMatrix;

uniform vec4 ambientMat;
uniform vec4 diffuseMat;
uniform vec4 specMat;

in vec3 position;
in vec3 vector;

out vec3 v;
out vec4 vertexColour;

void main()
{
    // Compute vector
    v = mat3(NormalMatrix) * vector;
    
    // Calculate colour
    vertexColour = ambientMat + diffuseMat + specMat;
    gl_Position = MVMatrix * vec4(position.xyz, 1.);
}
