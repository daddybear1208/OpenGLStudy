#include <iostream>
#ifdef _WIN32
#include <Windows.h>
#endif
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/GLU.h>
#include <GL/freeglut.h>

using namespace std;

static GLfloat vertex[] = { 10, 10, 0,
					200, 10, 0,
					200, 200, 0,
					10, 200, 0};

static GLfloat colors[] = {255,0,0,
						255,0,0,
						0,255,0,
						0,255,0};

GLbyte indices[] = { 0, 1, 2, 3 };

void init()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glVertexPointer(3, GL_FLOAT, 0, vertex);
	glColorPointer(3, GL_FLOAT, 0, colors);
	glDrawArrays(GL_QUADS, 0, 4);
	glutSwapBuffers();
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