#version 150

uniform vec4 ambientMat;
uniform vec4 diffuseMat;
uniform vec4 specMat;

uniform mat4 MVMatrix;
uniform mat4 MVPMatrix;
uniform mat4 NormalMatrix;

uniform vec3 lightPos;

in vec3 position;
in vec3 vector;

out vec4 colourV;

void main()
{
    // Define material specs
    float specPow = 5.;
    
    // Compute vectors
    vec4 p = MVMatrix * vec4(position.xyz, 1.);
    vec3 N = normalize(mat3(NormalMatrix) * vector);
    vec3 L = normalize(lightPos - p.xyz);
    vec3 E = normalize(-p.xyz);
    vec3 R = normalize(reflect(-L,N));
    
    // Calculate colour
    vec4 ambient = ambientMat;
    vec4 diffuse = clamp( diffuseMat * max(dot(N,L), 0.0)  , 0.0, 1.0 ) ;
    vec4 spec = clamp ( specMat * pow(max(dot(R,E),0.0), 0.3*specPow) , 0.0, 1.0 );
    colourV = ambient + diffuse + spec;
    
	gl_Position =  MVPMatrix * vec4(position.xyz, 1.);
}
