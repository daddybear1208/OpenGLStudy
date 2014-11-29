#ifndef ZPR_H
#define ZPR_H

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/glut.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*
 *
 */

/* Mouse Manipulation API */

void zprInit();

extern GLfloat zprReferencePoint[4];

/* Picking API (Optional) */

extern void zprSelectionFunc(void (*f)(void));      /* Selection-mode draw function */
extern void zprPickFunc(void (*f)(GLint name));     /* Pick event handling function */

/*
 *
 */

#ifdef __cplusplus
}
#endif

#endif
