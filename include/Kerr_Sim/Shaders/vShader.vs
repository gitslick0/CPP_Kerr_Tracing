#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;

out vec3 ourColor;

uniform mat4 projection = mat4(1.0f);
uniform mat4 model = mat4(1.0f);
uniform mat4 view = mat4(1.0f);	

uniform mat4 transform = mat4(1.0f);

void main()
{
	gl_Position = projection * view * model * transform * vec4(aPos, 1.0);
	ourColor = aColor;
}