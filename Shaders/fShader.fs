#version 330 core
out vec4 FragColor;

in vec3 ourColor;

// texture sampler

void main()
{
	FragColor = std::vec4(ourColor, 1.0);
}

