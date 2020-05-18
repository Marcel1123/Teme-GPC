//#define _USE_MATH_DEFINES
//#define _CRT_SECURE_NO_WARNINGS
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <assert.h>
//#include <float.h>
//#include<vector>
//#include "glut.h"
//using namespace std;
//#define dim 600
//
//#define NRITER_MB 50
//#define MODMAX_MB 10000000
//
//#define RX_MB 0.001
//#define RY_MB 0.001
//
//
//
//class CComplex {
//public:
//	CComplex() : re(0.0), im(0.0) {}
//	CComplex(double re1, double im1) : re(re1 * 1.0), im(im1 * 1.0) {}
//	CComplex(const CComplex& c) : re(c.re), im(c.im) {}
//	~CComplex() {}
//
//	CComplex& operator=(const CComplex& c)
//	{
//		re = c.re;
//		im = c.im;
//		return *this;
//	}
//
//	double getRe() { return re; }
//	void setRe(double re1) { re = re1; }
//
//	double getIm() { return im; }
//	void setIm(double im1) { im = im1; }
//
//	double getModul() { return sqrt(re * re + im * im); }
//
//	int operator==(CComplex& c1)
//	{
//		return ((re == c1.re) && (im == c1.im));
//	}
//
//	CComplex pow2()
//	{
//		CComplex rez;
//		rez.re = powl(re * 1.0, 2) - powl(im * 1.0, 2);
//		rez.im = 2.0 * re * im;
//		return rez;
//	}
//
//	friend CComplex operator+(const CComplex& c1, const CComplex& c2);
//	friend CComplex operator*(CComplex& c1, CComplex& c2);
//
//	void print(FILE* f)
//	{
//		fprintf(f, "%.20f%+.20f i", re, im);
//	}
//
//private:
//	double re, im;
//};
//
//CComplex operator+(const CComplex& c1, const CComplex& c2)
//{
//	CComplex rez(c1.re + c2.re, c1.im + c2.im);
//	return rez;
//}
//
//CComplex operator*(CComplex& c1, CComplex& c2)
//{
//	CComplex rez(c1.re * c2.re - c1.im * c2.im,
//		c1.re * c2.im + c1.im * c2.re);
//	return rez;
//}
//
//class Multime {
//
//public:
//    Multime()
//    {
//        m.nriter = NRITER_MB;
//        m.modmax = MODMAX_MB;
//    }
//
//    Multime(CComplex& c)
//    {
//        m.c = c;
//        m.nriter = NRITER_MB;
//        m.modmax = MODMAX_MB;
//    }
//
//    ~Multime() {}
//
//    void setmodmax(double v) { assert(v <= MODMAX_MB); m.modmax = v; }
//
//    double getmodmax() { return m.modmax; }
//
//    void setnriter(int v) { assert(v <= NRITER_MB); m.nriter = v; }
//
//    int getnriter() { return m.nriter; }
//
//    int isIn(CComplex& c)
//    {
//        int rez = 0;
//        CComplex z0, z_n;
//
//        z0 = CComplex(0.0, 0.0);
//        for (int i = 1; i < m.nriter; i++)
//        {
//            z_n = z0 * z0 + c;
//            if (z_n == z0)
//            {
//               return -1;
//            }
//            else if (z_n.getModul() > 2.0)
//            {
//               return i;
//            }
//            z0 = z_n;
//        }
//        return 0;
//    }
//    
//
//    void display(double xmin, double ymin, double xmax, double ymax)
//    {
//        glPushMatrix();
//        glLoadIdentity();
//        glTranslated((xmin + xmax) * 2.0 / (xmin - xmax), (ymin + ymax) * 2.0 / (ymin - ymax), 0);
//        glScaled(2.0 / (xmax - xmin), 2.0 / (ymax - ymin), 1);
//        glBegin(GL_POINTS);
//        for (double x = xmin; x <= xmax; x += RX_MB)
//            for (double y = ymin; y <= ymax; y += RY_MB)
//            {
//                CComplex z(x, y);
//                int r = isIn(z);
//
//                if (r == 0)
//                {
//                    glColor3f(1.0, 0.1, 0.1);
//                    glVertex2d(x, y);
//                }
//                else if (r == -1)
//                {
//                    glColor3f(0.1, 0.1, 1.0);
//                    glVertex2d(x, y);
//
//                }
//                else {
//                    glColor3f(0.1, 1.0, 0.1);
//                    glVertex2d(x, y);
//                }
//            }
//        fprintf(stdout, "STOP\n");
//        glEnd();
//        glPopMatrix();
//    }
//private:
//    struct SDate {
//        CComplex c;
//        int nriter;
//        double modmax;
//    } m;
//};
//
/////############################################################################
/////############################### EXERCITIUL 4 ###############################
/////############################################################################
//#define dim 600
//
//unsigned char prevKey;
//int nivel = 0;
//
//
//class C2coord
//{
//public:
//    C2coord()
//    {
//        m.x = m.y = 0;
//    }
//
//    C2coord(double x, double y)
//    {
//        m.x = x;
//        m.y = y;
//    }
//
//    C2coord(const C2coord& p)
//    {
//        m.x = p.m.x;
//        m.y = p.m.y;
//    }
//
//    C2coord& operator=(C2coord& p)
//    {
//        m.x = p.m.x;
//        m.y = p.m.y;
//        return *this;
//    }
//
//    int operator==(C2coord& p)
//    {
//        return ((m.x == p.m.x) && (m.y == p.m.y));
//    }
//
//protected:
//    struct SDate
//    {
//        double x, y;
//    } m;
//};
//
//class CPunct : public C2coord
//{
//public:
//    CPunct() : C2coord(0.0, 0.0)
//    {}
//
//    CPunct(double x, double y) : C2coord(x, y)
//    {}
//
//    CPunct& operator=(const CPunct& p)
//    {
//        m.x = p.m.x;
//        m.y = p.m.y;
//        return *this;
//    }
//
//    void getxy(double& x, double& y)
//    {
//        x = m.x;
//        y = m.y;
//    }
//
//    double getX()
//    {
//        return m.x;
//    }
//
//    double getY()
//    {
//        return m.y;
//    }
//
//    int operator==(CPunct& p)
//    {
//        return ((m.x == p.m.x) && (m.y == p.m.y));
//    }
//
//    void marcheaza()
//    {
//        glBegin(GL_POINTS);
//        glVertex2d(m.x, m.y);
//        glEnd();
//    }
//
//    void print(FILE* fis)
//    {
//        fprintf(fis, "(%+f,%+f)", m.x, m.y);
//    }
//
//};
//
//class CVector : public C2coord
//{
//public:
//    CVector() : C2coord(0.0, 0.0)
//    {
//        normalizare();
//    }
//
//    CVector(double x, double y) : C2coord(x, y)
//    {
//        normalizare();
//    }
//
//    CVector& operator=(CVector& p)
//    {
//        m.x = p.m.x;
//        m.y = p.m.y;
//        return *this;
//    }
//
//    int operator==(CVector& p)
//    {
//        return ((m.x == p.m.x) && (m.y == p.m.y));
//    }
//
//    CPunct getDest(CPunct& orig, double lungime)
//    {
//        double x, y;
//        orig.getxy(x, y);
//        CPunct p(x + m.x * lungime, y + m.y * lungime);
//        return p;
//    }
//
//    void rotatie(double grade)
//    {
//        double x = m.x;
//        double y = m.y;
//        double t = 2 * (4.0 * atan(1.0)) * grade / 360.0;
//        m.x = x * cos(t) - y * sin(t);
//        m.y = x * sin(t) + y * cos(t);
//        normalizare();
//    }
//
//    void deseneaza(CPunct p, double lungime)
//    {
//        double x, y;
//        p.getxy(x, y);
//        glColor3f(1.0, 0.1, 0.1);
//        glBegin(GL_LINE_STRIP);
//        glVertex2d(x, y);
//        glVertex2d(x + m.x * lungime, y + m.y * lungime);
//        glEnd();
//    }
//
//    void print(FILE* fis)
//    {
//        fprintf(fis, "%+fi %+fj", C2coord::m.x, C2coord::m.y);
//    }
//
//private:
//    void normalizare()
//    {
//        double d = sqrt(C2coord::m.x * C2coord::m.x + C2coord::m.y * C2coord::m.y);
//        if (d != 0.0)
//        {
//            C2coord::m.x = C2coord::m.x * 1.0 / d;
//            C2coord::m.y = C2coord::m.y * 1.0 / d;
//        }
//    }
//};
//
//class CImagine2
//{
//public:
//    void arborePerron(double lungime,
//        int nivel,
//        double factordiviziune,
//        CPunct p,
//        CVector v)
//    {
//        assert(factordiviziune != 0);
//        CPunct p1, p2;
//        if (nivel == 0)
//        {
//        }
//        else
//        {
//            v.rotatie(-45);
//            v.deseneaza(p, lungime);
//            p1 = v.getDest(p, lungime);
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//
//            v.rotatie(90);
//            v.deseneaza(p, lungime);
//            p1 = v.getDest(p, lungime);
//            p2 = p1;
//
//            v.rotatie(15);
//            v.deseneaza(p1, lungime);
//            p1 = v.getDest(p1, lungime);
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//
//            p1 = p2;
//            v.rotatie(-60);
//            v.deseneaza(p1, lungime);
//            p1 = v.getDest(p1, lungime);
//            p2 = p1;
//
//            v.rotatie(-90);
//            v.deseneaza(p1, lungime * (factordiviziune + 0.1));
//            p1 = v.getDest(p1, lungime * (factordiviziune + 0.1));
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//
//            p1 = p2;
//            v.rotatie(110);
//            v.deseneaza(p1, lungime * (factordiviziune + 0.1));
//            p1 = v.getDest(p1, lungime * (factordiviziune + 0.1));
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//        }
//    }
//
//    void afisare(double lungime, int nivel)
//    {
//        CVector v(0.0, -1.0);
//        CPunct p(0.0, 1.0);
//
//        v.deseneaza(p, 0.25);
//        p = v.getDest(p, 0.25);
//        arborePerron(lungime, nivel, 0.4, p, v);
//    }
//};
//
//class CPatrat {
//public:
//    CPunct ld, lu, ru, rd;
//    double length;
//
//    CPatrat(
//        CPunct ld,
//        CPunct lu,
//        CPunct ru,
//        CPunct rd,
//        double length
//    ) {
//        this->lu = lu;
//        this->ld = ld;
//        this->rd = rd;
//        this->ru = ru;
//        this->length = length;
//    }
//
//    CPatrat patratMijloc()
//    {
//        return CPatrat(
//            /// new ld point
//            CPunct(ld.getX() + length / 3, ld.getY() + length / 3),
//            /// new lu point
//            CPunct(ld.getX() + length / 3, ld.getY() + (length / 3) * 2),
//            /// new ru point
//            CPunct(ld.getX() + (length / 3) * 2, ld.getY() + (length / 3) * 2),
//            /// new rd point
//            CPunct(ld.getX() + (length / 3) * 2, ld.getY() + length / 3),
//            /// new length
//            length / 3
//        );
//    }
//
//    CPatrat mutarePatrat(CPatrat object_, double distantaX, double distantaY)
//    {
//        return CPatrat(
//            /// new ld point
//            CPunct(object_.ld.getX() + distantaX, object_.ld.getY() + distantaY),
//            /// new lu point
//            CPunct(object_.lu.getX() + distantaX, object_.lu.getY() + distantaY),
//            /// new ru point
//            CPunct(object_.ru.getX() + distantaX, object_.ru.getY() + distantaY),
//            /// new rd point
//            CPunct(object_.rd.getX() + distantaX, object_.rd.getY() + distantaY),
//            /// new length
//            object_.length
//        );
//    }
//
//    void draw()
//    {
//        glColor3f(1.0, 0.1, 0.1);
//        glBegin(GL_LINE_LOOP);
//
//        /// draw ld
//        glVertex2d(ld.getX(), ld.getY());
//
//        /// draw lu
//        glVertex2d(lu.getX(), lu.getY());
//
//        /// draw ru
//        glVertex2d(ru.getX(), ru.getY());
//
//        /// draw rd
//        glVertex2d(rd.getX(), rd.getY());
//
//        glEnd();
//        glFlush();
//    }
//};
//
//class CImagine1 {
//public:
//    void drawImagine1(CPatrat object_, int nivel)
//    {
//        if (nivel == 0)
//        {
//            object_.patratMijloc().draw();
//        }
//        else
//        {
//            CPatrat patrat = object_.patratMijloc();
//            patrat.draw();
//
//            /// Desenarea celor 8 patrate
//            drawImagine1(patrat.mutarePatrat(patrat, -patrat.length, patrat.length), nivel - 1);
//            drawImagine1(patrat.mutarePatrat(patrat, 0, patrat.length), nivel - 1);
//
//            drawImagine1(patrat.mutarePatrat(patrat, patrat.length, patrat.length), nivel - 1);
//            drawImagine1(patrat.mutarePatrat(patrat, patrat.length, 0), nivel - 1);
//            
//            drawImagine1(patrat.mutarePatrat(patrat, patrat.length, -patrat.length), nivel - 1);
//            drawImagine1(patrat.mutarePatrat(patrat, 0, -patrat.length), nivel - 1);
//            
//            drawImagine1(patrat.mutarePatrat(patrat, -patrat.length, -patrat.length), nivel - 1);
//            drawImagine1(patrat.mutarePatrat(patrat, -patrat.length, 0), nivel - 1);
//        }
//    }
//
//    void draw(int nivel)
//    {
//        double initialX, initialY, length;
//        initialX = -0.5;
//        initialY = -0.5;
//        length = 1;
//
//        CPatrat patratInitial = CPatrat(
//            CPunct(initialX, initialY), // ld
//            CPunct(initialX, initialY + length), // lu
//            CPunct(initialX + length, initialY + length), // ru
//            CPunct(initialX + length, initialY), // rd
//            length
//        );
//
//        patratInitial.draw();
//        drawImagine1(patratInitial, nivel);
//    }
//};
//
//void Display2() 
//{
//    CImagine1 object_;
//    object_.draw(nivel);
//
//    char cl[3];
//    sprintf(cl, "%2d", nivel);
//    glRasterPos2d(-0.98, -0.98);
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, cl[0]);
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, cl[1]);
//
//    glRasterPos2d(-1.0, -0.9);
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'I');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'm');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'a');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'g');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'n');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '1');
//
//    nivel++;
//}
//
//void Display3() {
//    CImagine2 cap;
//
//    char c[3];
//    sprintf(c, "%2d", nivel);
//    glRasterPos2d(-0.98, -0.98);
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);
//
//    glRasterPos2d(-1.0, -0.9);
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'I');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'm');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'a');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'g');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'n');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
//    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '2');
//
//    glPushMatrix();
//    glLoadIdentity();
//    glScaled(0.4, 0.4, 1);
//    glTranslated(-0.5, 1, 0.0);
//    cap.afisare(1, nivel);
//    glPopMatrix();
//    nivel++;
//}
//
//void Display1() {
//    CComplex c(-0.12375, 0.056805);
//    Multime cjf(c);
//
//    glColor3f(1.0, 0.1, 0.1);
//    cjf.setnriter(15);
//    cjf.display(-2,-2,2,2);
//}
//
//void Init(void) {
//
//    glClearColor(1.0, 1.0, 1.0, 1.0);
//
//    glLineWidth(1);
//
//    glPointSize(3);
//
//    glPolygonMode(GL_FRONT, GL_LINE);
//}
//
//void Display(void)
//{
//    switch (prevKey)
//    {
//    case '0':
//        glClear(GL_COLOR_BUFFER_BIT);
//        nivel = 0;
//        fprintf(stderr, "nivel = %d\n", nivel);
//        break;
//    case '1':
//        glClear(GL_COLOR_BUFFER_BIT);
//        Display1();
//        break;
//    case '2':
//        glClear(GL_COLOR_BUFFER_BIT);
//        Display2();
//        break;
//    case '3':
//        glClear(GL_COLOR_BUFFER_BIT);
//        Display3();
//        break;
//    default:
//        break;
//    }
//
//    glFlush();
//}
//
//void Reshape(int w, int h)
//{
//    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
//}
//
//void KeyboardFunc(unsigned char key, int x, int y)
//{
//    prevKey = key;
//    if (key == 27) // escape
//        exit(0);
//    glutPostRedisplay();
//}
//
//void MouseFunc(int button, int state, int x, int y)
//{
//}
//
//int main(int argc, char** argv)
//{
//    glutInit(&argc, argv);
//
//    glutInitWindowSize(dim, dim);
//
//    glutInitWindowPosition(100, 100);
//
//    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
//
//    glutCreateWindow(argv[0]);
//
//    Init();
//
//    glutReshapeFunc(Reshape);
//
//    glutKeyboardFunc(KeyboardFunc);
//
//    glutMouseFunc(MouseFunc);
//
//    glutDisplayFunc(Display);
//
//    glutMainLoop();
//
//    return 0;
//}