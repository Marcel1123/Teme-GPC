
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <algorithm>
#include <iostream>
#include <list>
#include "glut.h"

// dimensiunea ferestrei in pixeli
#define dim 300
using namespace std;

unsigned char prevKey;

void Display4() {
    double xmax, ymax, xmin, ymin;
    double a = 0.3, b = 0.2;
    double pi = 4 * atan(1.0);
    double ratia = 0.01;
    
    xmax = -FLT_MAX;
    xmin = FLT_MAX;

    ymax = -FLT_MAX;
    ymin = FLT_MAX;

    for (double t = -pi + ratia; t < pi; t += ratia) {
        double x1, y1;

        x1 = 2 * (a * cos(t) + b) * cos(t);
        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        y1 = 2 * (a * cos(t) + b) * sin(t);
        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    // afisarea punctelor propriu-zise precedata de scalare
    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (double t = -pi + ratia; t < pi; t += ratia) {
        double x1, y1;

        x1 = (2 * (a * cos(t) + b) * cos(t)) / (xmax * 1.2);
        y1 = (2 * (a * cos(t) + b) * sin(t)) / (ymax * 1.8);

        glVertex2f(x1, y1);
    }
    glEnd();
}

void Display1() {
   double xmax, ymax, xmin, ymin;
   double a = 1, b = 2;
   double pi = 4 * atan(1.0);
   double ratia = 0.05;
   double t;

   // calculul valorilor maxime/minime ptr. x si y
   // aceste valori vor fi folosite ulterior la scalare
   xmax = a - b - 1;
   xmin = a + b + 1;
   ymax = ymin = 0;
   for (t = - pi/2 + ratia; t < pi / 2; t += ratia) {
      double x1, y1, x2, y2;
      x1 = a + b * cos(t);
      xmax = (xmax < x1) ? x1 : xmax;
      xmin = (xmin > x1) ? x1 : xmin;

      x2 = a - b * cos(t);
      xmax = (xmax < x2) ? x2 : xmax;
      xmin = (xmin > x2) ? x2 : xmin;

      y1 = a * tan(t) + b * sin(t);
      ymax = (ymax < y1) ? y1 : ymax;
      ymin = (ymin > y1) ? y1 : ymin;

      y2 = a * tan(t) - b * sin(t);
      ymax = (ymax < y2) ? y2 : ymax;
      ymin = (ymin > y2) ? y2 : ymin;
   }

   xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
   ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

   // afisarea punctelor propriu-zise precedata de scalare
   glColor3f(1,0.1,0.1); // rosu
   glBegin(GL_LINE_STRIP); 
   for (t = - pi/2 + ratia; t < pi / 2; t += ratia) {
      double x1, y1, x2, y2;
      x1 = (a + b * cos(t)) / xmax;
      x2 = (a - b * cos(t)) / xmax;
      y1 = (a * tan(t) + b * sin(t)) / ymax;
      y2 = (a * tan(t) - b * sin(t)) / ymax;

      glVertex2f(x1,y1);
   }
   glEnd();

   glBegin(GL_LINE_STRIP); 
   for (t = - pi/2 + ratia; t < pi / 2; t += ratia) {
      double x1, y1, x2, y2;
      x1 = (a + b * cos(t)) / xmax;
      y1 = (a * tan(t) + b * sin(t)) / ymax;
      x2 = (a - b * cos(t)) / xmax;
      y2 = (a * tan(t) - b * sin(t)) / ymax;

      glVertex2f(x2,y2);
   }
   glEnd();
}

void Display2() {
   double pi = 4 * atan(1.0);
   double ratia = 0.05;
   double xmax = 8 * pi;
   double ymax = exp(1.1);

   // afisarea punctelor propriu-zise precedata de scalare
   glColor3f(1,0.1,0.1); // rosu
   glBegin(GL_LINE_STRIP); 
   for (double x = 0; x < xmax; x += ratia) {
      double x1, y1;
      x1 = x / xmax;
      y1 = (fabs(sin(x)) * exp(-sin(x))) / ymax;

      glVertex2f(x1,y1);
   }
   glEnd();
}

void Display3() {
    double xmax = 100;
    double xmin = 0;
    double ymax = -FLT_MAX;
    double ymin = FLT_MAX; 
    double ratia = 0.05;

    for (double i = 0; i <= xmax; i += ratia) {
        double local_y;
        double left_nuber = floor(i);
        double right_number = ceil(i);
        if (i == 0) {
            local_y = 1;
        }
        else if (i - left_nuber < right_number - i) {
            local_y = (i - left_nuber) / i;
        }
        else {
            local_y = (right_number - i) / i;
        }
        ymax = (ymax < local_y) ? local_y : ymax;
        ymin = (ymin > local_y) ? local_y : ymin;
    }
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    glVertex2d(0, 1 / (ymax * 1.6));
    for (double i = 0 + ratia; i <= xmax; i += ratia) {
        double left_nuber = floor(i);
        double right_number = ceil(i);
        double local_y = 0;

        if (i == 0) {
            local_y = 1;
        }
        else if (i - left_nuber < right_number - i) {
            local_y = (i - left_nuber) / i;
        }
        else {
            local_y = (right_number - i) / i;
        }
        glVertex2f(i/ (xmax * 1.6), local_y / (ymax * 1.6));
    }
    glEnd();
}

void Display6() {
    double xmax, ymax, xmin, ymin;
    double a = 0.1, b = 0.2;
    double pi = 4 * atan(1.0);
    double ratia = 0.05;

    xmax = -FLT_MAX;
    xmin = FLT_MAX;

    ymax = -FLT_MAX;
    ymin = FLT_MAX;

    for (double t = -2 * pi + ratia; t < 2 * pi; t += ratia) {
        double x1, y1;
        x1 = a * t - b * sin(t);
        y1 = a - b * cos(t);

        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;

    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (double t = -3 * pi + ratia; t < 3*pi; t += ratia) {
        double x1, y1;

        x1 = (a * t - b * sin(t)) / (xmax* 1.4);
        y1 = (a - b * cos(t)) / (ymax * 3);

        glVertex2f(x1, y1);
    }
    glEnd();
}

void Init(void) {

   glClearColor(1.0,1.0,1.0,1.0);

   glLineWidth(1);

//   glPointSize(4);

   glPolygonMode(GL_FRONT, GL_LINE);
}

void Display7() {
    double R = 0.1;
    double r = 0.3;
    double xmax, ymax, xmin, ymin;
    double pi = 4 * atan(1.0);
    double ratia = 0.05;

    xmax = -FLT_MAX;
    xmin = FLT_MAX;

    ymax = -FLT_MAX;
    ymin = FLT_MAX;

    for (double t = 0; t <= 2 * pi; t += ratia) {
        double x1, y1;
        x1 = (R + r) * cos((r / R) * t) - r * cos(t + (r / R) * t);
        y1 = (R + r) * sin((r / R) * t) - r * sin(t + (r / R) * t);

        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (double t = 0; t <= 2 * pi; t += ratia) {
        double x1, y1;

        x1 = (R + r) * cos((r / R) * t) - r * cos(t + (r / R) * t);
        y1 = (R + r) * sin((r / R) * t) - r * sin(t + (r / R) * t);

        glVertex2f(x1 / (xmax * 1.3), y1 / (ymax * 1.3));
    }
    glEnd();
}

void Display8() {
    double R = 0.1;
    double r = 0.3;
    double xmax, ymax, xmin, ymin;
    double pi = 4 * atan(1.0);
    double ratia = 0.05;

    xmax = -FLT_MAX;
    xmin = FLT_MAX;

    ymax = -FLT_MAX;
    ymin = FLT_MAX;

    for (double t = 0; t <= 2 * pi; t += ratia) {
        double x1, y1;
        x1 = (R - r) * cos((r / R) * t) - r * cos(t - (r / R) * t);
        y1 = (R - r) * sin((r / R) * t) - r * sin(t - (r / R) * t);

        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    glColor3f(1, 0.1, 0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (double t = 0; t <= 2 * pi; t += ratia) {
        double x1, y1;

        x1 = (R - r) * cos((r / R) * t) - r * cos(t - (r / R) * t);
        y1 = (R - r) * sin((r / R) * t) - r * sin(t - (r / R) * t);

        glVertex2f(x1 / (xmax * 1.6), y1 / (ymax * 1.6));
    }
    glEnd();
}

void Display9() {
    double pi = 4 * atan(1.0);
    double ratia = 0.01;
    double xmax = -FLT_MAX;
    double ymax = -FLT_MAX;
    double xmin = FLT_MAX;
    double ymin = FLT_MAX;
    double a = 0.4;

    for (double t = -(pi / 4) + ratia; t < pi / 4; t += ratia) {
        double r1, r2, x1, y1, x2, y2;

        r1 = -a * (sqrt(2 * cos(2 * t)));
        x1 = r1 * cos(t);
        y1 = r1 * sin(t);

        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;

        r2 = a * (sqrt(2 * cos(2 * t)));
        x2 = r2 * cos(t);
        y2 = r2 * sin(t);

        xmax = (xmax < x2) ? x2 : xmax;
        xmin = (xmin > x2) ? x2 : xmin;

        ymax = (ymax < y2) ? y2 : ymax;
        ymin = (ymin > y2) ? y2 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    glColor3f(1, 0.1, 0.1);
    glBegin(GL_LINE_STRIP);
    for (double t = (pi / 4) - (ratia / 24); t >= -pi / 4; t -= ratia) {
        double r1, x1, y1;

        r1 = - a * (sqrt(2 * cos(2 * t)));
        x1 = r1 * cos(t);
        y1 = r1 * sin(t);

        glVertex2f(x1 / (xmax * 1.3), y1 / (ymax * 3.3));
    }

    for (double t = -(pi / 4); t < (pi / 4); t += ratia) {
        double r2, x2, y2;

        r2 =  a * (sqrt(2 * cos(2 * t)));
        x2 = r2 * cos(t);
        y2 = r2 * sin(t);

        glVertex2f(x2 / (xmax * 1.3), y2 / (ymax * 3.3));
    }
    glEnd();
}

void Display10() {
    double pi = 4 * atan(1.0);
    double ratia = 0.05;
    double xmax = -FLT_MAX;
    double ymax = -FLT_MAX;
    double xmin = FLT_MAX;
    double ymin = FLT_MAX;
    double a = 0.02;

    for (double t = ratia; t < pi; t += ratia) {
        double r1, x1, y1;

        r1 = a * exp(1 + t);
        x1 = r1 * cos(t);
        y1 = r1 * sin(t);

        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    glColor3f(1, 0.1, 0.1);
    glBegin(GL_LINE_STRIP);

    for (double t = ratia; t < pi; t += ratia) {
        double r1, x1, y1;

        r1 = a * exp(1 + t);
        x1 = r1 * cos(t);
        y1 = r1 * sin(t);
        glVertex2f(x1 / (xmax * 1.1), y1 / (ymax * 3.1));
    }
    glEnd();
}

void Display5() {
    double pi = 4 * atan(1.0);
    double ratia = 0.05;
    double xmax = -FLT_MAX;
    double ymax = -FLT_MAX;
    double xmin = FLT_MAX;
    double ymin = FLT_MAX;
    double a = 0.2;

    for (double t = -(pi / 2) + ratia; t < -(pi / 6); t += ratia)
    {
        if (t != -(pi / 6) || t != (pi / 6)) {
            double x1, y1;

            x1 = a / (4 * pow(cos(t), 2) - 3);
            y1 = (a * tan(t)) / (4 * pow(cos(t), 2) - 3);

            xmax = (xmax < x1) ? x1 : xmax;
            xmin = (xmin > x1) ? x1 : xmin;

            ymax = (ymax < y1) ? y1 : ymax;
            ymin = (ymin > y1) ? y1 : ymin;
        }
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    glColor3f(1, 0.1, 0.1);
    glBegin(GL_LINE_STRIP);
    bool start = true;
    double aux = 0.0, auy = 0.0;
    int position = 0;

    for (double t = -(pi / 2) + ratia; t < -(pi / 6); t += ratia) {
        if (t != -(pi / 6) || t != (pi / 6)) {
            double x1, y1;

            x1 = a / (4 * pow(cos(t), 2) - 3);
            y1 = (a * tan(t)) / (4 * pow(cos(t), 2) - 3);

            if (start == true) {
                aux = x1 / (xmax * 1.2);
                auy = y1 / (ymax * 1.2);
                start = false;
            }

            glVertex2f(x1 / (xmax * 1.2), y1 / (ymax * 1.2));
        }
    }
    glVertex2f(-0.8, 0.827);
    glVertex2f(aux, auy);
    glEnd();
    
    glColor3f(1, 0.1, 0.1);
    glBegin(GL_TRIANGLES);
    {
        glVertex2f(-0.8, 0.827);
        glVertex2f(aux, auy - 0.06);
        glVertex2f(aux, auy - 0.16);

        glVertex2f(-0.8, 0.827);
        glVertex2f(aux, auy - 0.21);
        glVertex2f(aux, auy - 0.29);

        glVertex2f(-0.8, 0.827);
        glVertex2f(aux, auy - 0.33);
        glVertex2f(aux, auy - 0.4);

        glVertex2f(-0.8, 0.827);
        glVertex2f(aux, auy - 0.44);
        glVertex2f(aux, auy - 0.5);

        glVertex2f(-0.8, 0.827);
        glVertex2f(aux - 0.001, auy - 0.54);
        glVertex2f(aux - 0.001, auy - 0.58);
    }
    glEnd();
}

void Display(void) {
   glClear(GL_COLOR_BUFFER_BIT);

   switch(prevKey) {
   case '1':
      Display1();
      break;
   case '2':
      Display2();
      break;
   case '3':
       Display3();
       break;
   case '4':
       Display4();
       break;
   case '5':
       Display5();
       break;
   case '6':
       Display6();
       break;
   case '7':
       Display7();
       break;
   case '8':
       Display8();
       break;
   case '9':
       Display9();
       break;
   case '0':
       Display10();
       break;
   default:
      break;
   }

   glFlush();
}

void Reshape(int w, int h) {
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

void KeyboardFunc(unsigned char key, int x, int y) {
   prevKey = key;
   if (key == 27) // escape
      exit(0);
   glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
}

int main(int argc, char** argv) {

   glutInit(&argc, argv);
   
   glutInitWindowSize(dim, dim);

   glutInitWindowPosition(100, 100);

   glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);

   glutCreateWindow (argv[0]);

   Init();

   glutReshapeFunc(Reshape);
   
   glutKeyboardFunc(KeyboardFunc);
   
   glutMouseFunc(MouseFunc);

   glutDisplayFunc(Display);
   
   glutMainLoop();

   return 0;
}
