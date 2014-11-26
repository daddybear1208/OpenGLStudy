#include <iostream>
#ifdef _WIN32
#include <Windows.h>
#endif
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GL/freeglut.h>

using namespace std;

enum Buffer_IDs
{
	VerticesBuf,
	IndicesBuf,
	NumBuffers
};

GLuint buffers[NumBuffers];

GLfloat vertices[][3] = 
{
	{-100.0,-100.0,-100.0},
	{100.0,-100.0,-100.0},
	{100.0,100.0,-100.0},
	{-100.0,100.0,-100.0},
	{-100.0,-100.0,100.0},
	{100.0,-100.0,100.0},
	{100.0,100.0,100.0},
	{-100.0,100.0,100.0}
};

GLubyte indices[][4] = 
{
	{0,1,2,3},
	{4,7,6,5},
	{0,4,5,1},
	{3,2,6,7},
	{0,3,7,4},
	{1,5,6,2}
};

void init()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnableClientState(GL_VERTEX_ARRAY);
	glGenBuffers(NumBuffers, buffers);
	glBindBuffer(GL_ARRAY_BUFFER, buffers[VerticesBuf]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[IndicesBuf]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), 0, GL_STATIC_DRAW);
}

void display()
{
	glColor3f(1.0, 0.0, 0.0);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glDrawElements(GL_QUADS, 24, GL_UNSIGNED_BYTE, 0);
}

void reshape(int w, int h)
{
	glViewport(0, 0, GLsizei(w), GLsizei(h));
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, GLdouble(w), 0.0, GLdouble(h));
}


int main(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(512, 512);
	glutInitWindowPosition(100, 100);
	//glutInitContextVersion(4, 3);
	//glutInitContextProfile(GLUT_CORE_PROFILE);
	glutCreateWindow("OpenGL");
	if (glewInit())
	{
		std::cerr << "Unable to initialize GLEW, exiting" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMainLoop();
	return 0;
}