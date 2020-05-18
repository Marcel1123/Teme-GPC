#include <iostream>
#include "glut.h"

using namespace std;

unsigned char prevKey;

class GrilaCarteziana {

	int Numar_linii;
	int Numar_coloane;
	float epsilon = 0.3;
	float dl;
	float dc;

public:

	GrilaCarteziana(int linii, int coloane) {
		SetNumarLinii(linii);
		SetNumarColoane(coloane);
	}
	void SetNumarLinii(int linii) {
		this->Numar_linii = linii;
	}
	int GetNumarLinii() {
		return this->Numar_linii;
	}
	void SetNumarColoane(int coloane) {
		this->Numar_coloane = coloane;
	}
	int GetNumarColoane() {
		return this->Numar_coloane;
	}
	void DisplayGrila() {
		dc = (2 - 2 * epsilon) / (Numar_coloane - 1);
		dl = (2 - 2 * epsilon) / (Numar_linii - 1);
		int i = 0;

		
		for (i = 0; i < Numar_linii; i++) {
			glColor3f(0, 0, 0);
			glBegin(GL_LINES);
			glVertex2d(-1 + epsilon, -1  +  epsilon + i*dl );
			glVertex2d(1 - epsilon, -1 + epsilon + i*dl);
			glEnd();
		}
		
		for (i = 0; i < Numar_coloane; i += 1) {
			glColor3f(0, 0, 0);
			glBegin(GL_LINES);
			glVertex2d(-1 + epsilon + i*dc, -1 + epsilon);
			glVertex2d(-1 + epsilon + i*dc, 1 - epsilon);
			glEnd();
		}

		drawLine(0, 15, 15, 8, 1);
		drawLine(0, 0, 15, 5, 3);

	}

	void writePixel(int linie, int coloana) {
		GLfloat x =  -1 + epsilon + coloana * dc;
		GLfloat y = 1 - epsilon - linie * dl;
		drawFilledCircle(x, y, 0.03);
	}

	void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius) {
		int i;
		int triangleAmount = 50;

		GLfloat twicePi = 2.0f * 3.14;

		glColor3f(0.65, 0.65, 0.65);
		glBegin(GL_TRIANGLE_FAN);
		glVertex2f(x, y); 
		for (i = 0; i <= triangleAmount; i++) {
			glVertex2f(
				x + (radius * cos(i * twicePi / triangleAmount)),
				y + (radius * sin(i * twicePi / triangleAmount))
			);
		}
		glEnd();
	}

	void drawLine(int x0, int y0, int xn, int yn, int grosime) {
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		glVertex2d(-1 + epsilon + x0 * dc, 1 - epsilon - y0*dl);
		glVertex2d(-1 + epsilon + xn * dc, 1 - epsilon - yn*dl);
		glEnd();
		int pixeli_extra = grosime / 2;
		float y = 0;
		for (int x = x0; x <= xn; x += 1) {
			y = (float) (x - x0) * (yn - y0) / (xn - x0) + y0;
			writePixel(round(y), x);// 3.8 + 3.7+2.75+3 
			for (int j = 1; j <= pixeli_extra; j++) {
				if (round(y) + j <= Numar_linii) {
					writePixel(round(y) + j, x);
				}
				if (round(y) - j >= 0) {
					writePixel(round(y) - j, x);
				}
			}

		}
		
	}

};

GrilaCarteziana grila(16, 16);


void Init(void) {

	glClearColor(1.0, 1.0, 1.0, 1.0);

	// grosimea liniilor
	glLineWidth(3);

	// dimensiunea punctelor
	glPointSize(4);

	glPolygonMode(GL_FRONT, GL_LINE);
}

void Display(void) {
	printf("Call Display\n");

	// sterge buffer-ul indicat
	glClear(GL_COLOR_BUFFER_BIT);

	switch (prevKey) {
	case '1':
		grila.DisplayGrila();
		break;
	default:
		break;
	}

	glFlush();
}

void Reshape(int w, int h) {
	printf("Call Reshape : latime = %d, inaltime = %d\n", w, h);
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
}

void KeyboardFunc(unsigned char key, int x, int y) {
	printf("Ati tastat <%c>. Mouse-ul este in pozitia %d, %d.\n",
		key, x, y);
	prevKey = key;
	if (key == 27) // escape
		exit(0);
	glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
	printf("Call MouseFunc : ati %s butonul %s in pozitia %d %d\n",
		(state == GLUT_DOWN) ? "apasat" : "eliberat",
		(button == GLUT_LEFT_BUTTON) ?
		"stang" :
		((button == GLUT_RIGHT_BUTTON) ? "drept" : "mijlociu"),
		x, y);
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(100, 100);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutCreateWindow(argv[0]);
	Init();
	glutReshapeFunc(Reshape);
	glutKeyboardFunc(KeyboardFunc);
	glutMouseFunc(MouseFunc);
	glutDisplayFunc(Display);
	glutMainLoop();

	return 0;
}




