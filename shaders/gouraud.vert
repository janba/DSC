#version 150

uniform mat4 MVMatrix;
uniform mat4 MVPMatrix;
uniform mat4 NormalMatrix;
uniform vec3 lightPos;
uniform vec3 eyePos;

in vec3 position;
in vec3 normal;

out vec4 colourV;

void main()
{
    // Define material specs
    float specPow = 5.;
    vec4 ambientMat = vec4(0.3, 0.3, 0.3, 1.);
    vec4 diffuseMat = vec4(0.5, 0.5, .5, 1.);
    vec4 specMat = vec4(0.2, 0.2, 0.2, 1.);
    
    // Compute vectors
    vec4 p = MVMatrix * vec4(position.xyz, 1.);
    vec3 N = normalize(mat3(NormalMatrix) * normal);
    vec3 L = normalize(lightPos - p.xyz);
    vec3 E = normalize(eyePos - p.xyz);
    vec3 R = normalize(reflect(-L,N));
    
    // Calculate colour
    vec4 ambient = ambientMat;
    vec4 diffuse = clamp( diffuseMat * max(dot(N,L), 0.0)  , 0.0, 1.0 ) ;
    vec4 spec = clamp ( specMat * pow(max(dot(R,E),0.0),0.3*specPow) , 0.0, 1.0 );
    colourV = ambient + diffuse + spec;
    
    // Calculate position
    gl_Position = MVPMatrix * vec4(position.xyz, 1.);
}
