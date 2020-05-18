#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <vector>
#include "glut.h"

// dimensiunea ferestrei in pixeli
//#define dim 300
//#define deg_60 M_PI / 3
//#define deg_90 M_PI / 2
//#define deg_30 M_PI / 6
//using namespace std;
//unsigned char prevKey;
//int nivel = 0;
//
//// numarul maxim de iteratii pentru testarea apartenentei la mult.Julia-Fatou
//#define NRITER_JF 5000
//#define NRITER_M 5000
//// modulul maxim pentru testarea apartenentei la mult.Julia-Fatou
//#define MODMAX_JF 10000000
//// ratii ptr. CJuliaFatou
//#define RX_JF 0.01
//#define RY_JF 0.01
//#define RX_M 0.01
//#define RY_M 0.01
//
//class CComplex {
//public:
//    CComplex() : re(0.0), im(0.0) {}
//    CComplex(double re1, double im1) : re(re1 * 1.0), im(im1 * 1.0) {}
//    CComplex(const CComplex& c) : re(c.re), im(c.im) {}
//    ~CComplex() {}
//
//    CComplex& operator=(const CComplex& c)
//    {
//        re = c.re;
//        im = c.im;
//        return *this;
//    }
//
//    double getRe() { return re; }
//    void setRe(double re1) { re = re1; }
//
//    double getIm() { return im; }
//    void setIm(double im1) { im = im1; }
//
//    double getModul() { return sqrt(re * re + im * im); }
//
//    int operator==(CComplex& c1)
//    {
//        return ((re == c1.re) && (im == c1.im));
//    }
//
//    CComplex pow2()
//    {
//        CComplex rez;
//        rez.re = powl(re * 1.0, 2) - powl(im * 1.0, 2);
//        rez.im = 2.0 * re * im;
//        return rez;
//    }
//
//    friend CComplex operator+(const CComplex& c1, const CComplex& c2);
//    friend CComplex operator*(CComplex& c1, CComplex& c2);
//
//    void print(FILE* f)
//    {
//        fprintf(f, "%.20f%+.20f i", re, im);
//    }
//
//private:
//    double re, im;
//};
//
//CComplex operator+(const CComplex& c1, const CComplex& c2)
//{
//    CComplex rez(c1.re + c2.re, c1.im + c2.im);
//    return rez;
//}
//
//CComplex operator*(CComplex& c1, CComplex& c2)
//{
//    CComplex rez(c1.re * c2.re - c1.im * c2.im,
//        c1.re * c2.im + c1.im * c2.re);
//    return rez;
//}
//
//class CJuliaFatou {
//public:
//    CJuliaFatou()
//    {
//        // m.c se initializeaza implicit cu 0+0i
//
//        m.nriter = NRITER_JF;
//        m.modmax = MODMAX_JF;
//    }
//
//    CJuliaFatou(CComplex& c)
//    {
//        m.c = c;
//        m.nriter = NRITER_JF;
//        m.modmax = MODMAX_JF;
//    }
//
//    ~CJuliaFatou() {}
//
//    void setmodmax(double v) { assert(v <= MODMAX_JF); m.modmax = v; }
//    double getmodmax() { return m.modmax; }
//
//    void setnriter(int v) { assert(v <= NRITER_JF); m.nriter = v; }
//    int getnriter() { return m.nriter; }
//
//    // testeaza daca x apartine multimii Julia-Fatou Jc
//    // returneaza 0 daca apartine, -1 daca converge finit, +1 daca converge infinit
//    int isIn(CComplex& x)
//    {
//        int rez = 0;
//        // tablou in care vor fi memorate valorile procesului iterativ z_n+1 = z_n * z_n + c
//        CComplex z0, z1;
//
//        z0 = x;
//        for (int i = 1; i < m.nriter; i++)
//        {
//            z1 = z0 * z0 + m.c;
//            if (z1 == z0)
//            {
//                // x nu apartine m.J-F deoarece procesul iterativ converge finit
//                rez = -1;
//                break;
//            }
//            else if (z1.getModul() > m.modmax)
//            {
//                // x nu apartine m.J-F deoarece procesul iterativ converge la infinit
//                rez = 1;
//                break;
//            }
//            z0 = z1;
//        }
//
//        return rez;
//    }
//
//    // afisarea multimii J-F care intersecteaza multimea argument
//    void display(double xmin, double ymin, double xmax, double ymax)
//    {
//        glPushMatrix();
//        glLoadIdentity();
//
//        //   glTranslated((xmin + xmax) * 1.0 / (xmin - xmax), (ymin + ymax)  * 1.0 / (ymin - ymax), 0);
//        //    glScaled(1.0 / (xmax - xmin), 1.0 / (ymax - ymin), 1);
//            // afisarea propriu-zisa
//        glBegin(GL_POINTS);
//        for (double x = xmin; x <= xmax; x += RX_JF)
//            for (double y = ymin; y <= ymax; y += RY_JF)
//            {
//                CComplex z(x, y);
//                int r = isIn(z);
//                //        z.print(stdout);
//                if (r == 0)
//                {
//                    //          fprintf(stdout, "   \n");
//                    glVertex3d(x, y, 0);
//                }
//                else if (r == -1)
//                {
//                    //          fprintf(stdout, "   converge finit\n");
//                }
//                else if (r == 1)
//                {
//                    //          fprintf(stdout, "   converge infinit\n");
//                }
//            }
//        fprintf(stdout, "STOP\n");
//        glEnd();
//
//        glPopMatrix();
//    }
//
//private:
//    struct SDate {
//        CComplex c;
//        // nr. de iteratii
//        int nriter;
//        // modulul maxim
//        double modmax;
//    } m;
//};
//
//class Mandelbrot {
//private:
//    int nriter;
//public:
//
//    Mandelbrot()
//    {
//        nriter = NRITER_M;
//    }
//
//    ~Mandelbrot() {}
//
//    void setNrIter(int nriter)
//    {
//        assert(nriter <= NRITER_M);
//        this->nriter = nriter;
//    }
//
//    int getNrIter()
//    {
//        return nriter;
//    }
//
//    int isIn(CComplex& c)
//    {
//        CComplex z0 = CComplex(0, 0);
//
//        for (int i = 1; i <= nriter; i++)
//        {
//            CComplex z1 = z0.pow2() + c;
//
//            if (z1.getModul() > 2)
//            {
//                return i;
//            }
//
//            z0 = z1;
//        }
//
//        return 0;
//    }
//
//    void display(double xmin, double ymin, double xmax, double ymax)
//    {
//        glPushMatrix();
//        glLoadIdentity();
//        glBegin(GL_POINTS);
//        //glTranslated((xmin + xmax) * 1.0 / (xmin - xmax), (ymin + ymax)  * 1.0 / (ymin - ymax), 0);
//        //glScaled(1.0 / (xmax - xmin), 1.0 / (ymax - ymin), 1);
//
//        for (double x = xmin; x <= xmax; x += RX_M)
//            for (double y = ymin; y <= ymax; y += RY_M)
//            {
//                CComplex c = CComplex(x, y);
//                int result = isIn(c);
//
//                if (result == 0)
//                {
//                    glVertex3d((x + 0.5) / 2, y / 2, 0);
//                }
//                else if (result > 0)
//                {
//                    glColor3f(1.0 - (1.0 / result), 0.1 + (0.1 * result), 0.1);
//
//                    glVertex3d((x + 0.5) / 2, y / 2, 0);
//
//                    glColor3f(1.0, 0.1, 0.1);
//                }
//            }
//        fprintf(stdout, "STOP\n");
//        glEnd();
//
//        glPopMatrix();
//    }
//
//};
//
//// multimea Julia-Fatou pentru z0 = 0 si c = -0.12375+0.056805i
//void Display8() {
//    CComplex c(-0.12375, 0.056805);
//    CJuliaFatou cjf(c);
//
//    glColor3f(1.0, 0.1, 0.1);
//    cjf.setnriter(30);
//    cjf.display(-0.8, -0.4, 0.8, 0.4);
//}
//
//// multimea Julia-Fatou pentru z0 = 0 si c = -0.012+0.74i
//void Display9() {
//    CComplex c(-0.012, 0.74);
//    CJuliaFatou cjf(c);
//
//    glColor3f(1.0, 0.1, 0.1);
//    cjf.setnriter(30);
//    cjf.display(-1, -1, 1, 1);
//}
//
//void DisplayMandelbrot()
//{
//    Mandelbrot m = Mandelbrot();
//
//    glColor3f(1.0, 0.1, 0.1);
//    m.setNrIter(30);
//    m.display(-2, -2, 2, 2);
//}
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
//class CCurbaKoch
//{
//public:
//    void segmentKoch(double lungime, int nivel, CPunct& p, CVector v)
//    {
//        CPunct p1;
//        if (nivel == 0)
//        {
//            v.deseneaza(p, lungime);
//        }
//        else
//        {
//            //    v.print(stderr);
//            //    fprintf(stderr, "\n");
//            segmentKoch(lungime / 3.0, nivel - 1, p, v);
//            p1 = v.getDest(p, lungime / 3.0);
//            v.rotatie(60);
//            //    v.print(stderr);
//            //    fprintf(stderr, "\n");
//            segmentKoch(lungime / 3.0, nivel - 1, p1, v);
//            p1 = v.getDest(p1, lungime / 3.0);
//            v.rotatie(-120);
//            //    v.print(stderr);
//            //    fprintf(stderr, "\n");
//            segmentKoch(lungime / 3.0, nivel - 1, p1, v);
//            p1 = v.getDest(p1, lungime / 3.0);
//            v.rotatie(60);
//            //    v.print(stderr);
//            //    fprintf(stderr, "\n");
//            segmentKoch(lungime / 3.0, nivel - 1, p1, v);
//        }
//    }
//
//    void afisare(double lungime, int nivel)
//    {
//        CVector v1(sqrt(3.0) / 2.0, 0.5);
//        CPunct p1(-1.0, 0.0);
//
//        CVector v2(0.0, -1.0);
//        CPunct p2(0.5, sqrt(3.0) / 2.0);
//
//        CVector v3(-sqrt(3.0) / 2.0, 0.5);
//        CPunct p3(0.5, -sqrt(3.0) / 2.0);
//
//        segmentKoch(lungime, nivel, p1, v1);
//        segmentKoch(lungime, nivel, p2, v2);
//        segmentKoch(lungime, nivel, p3, v3);
//    }
//};
//
//class CArboreBinar
//{
//public:
//    void arboreBinar(double lungime, int nivel, CPunct& p, CVector v)
//    {
//        CPunct p1;
//        if (nivel == 0)
//        {
//            v.deseneaza(p, lungime);
//        }
//        else
//        {
//            arboreBinar(lungime, nivel - 1, p, v);
//            p1 = v.getDest(p, lungime);
//
//            v.rotatie(-45);
//            arboreBinar(lungime / 2.0, nivel - 1, p1, v);
//
//            v.rotatie(90);
//            arboreBinar(lungime / 2.0, nivel - 1, p1, v);
//        }
//    }
//
//    void afisare(double lungime, int nivel)
//    {
//        CVector v(0.0, -1.0);
//        CPunct p(0.0, 1.0);
//
//        arboreBinar(lungime, nivel, p, v);
//    }
//};
//
//class CArborePerron
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
//            v.rotatie(30);
//            v.deseneaza(p, lungime);
//            p1 = v.getDest(p, lungime);
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//
//            v.rotatie(-90);
//            v.deseneaza(p, lungime);
//            p1 = v.getDest(p, lungime);
//            p2 = p1;
//
//            v.rotatie(-30);
//            v.deseneaza(p1, lungime);
//            p1 = v.getDest(p1, lungime);
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//
//            p1 = p2;
//            v.rotatie(90);
//            v.deseneaza(p1, lungime);
//            p1 = v.getDest(p1, lungime);
//            p2 = p1;
//
//            v.rotatie(30);
//            v.deseneaza(p1, lungime);
//            p1 = v.getDest(p1, lungime);
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//
//            p1 = p2;
//            v.rotatie(-90);
//            v.deseneaza(p1, lungime);
//            p1 = v.getDest(p1, lungime);
//            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);
//        }
//    }
//
//    void afisare(double lungime, int nivel)
//    {
//        CVector v(0.0, 1.0);
//        CPunct p(0.0, -1.0);
//
//        v.deseneaza(p, 0.25);
//        p = v.getDest(p, 0.25);
//        arborePerron(lungime, nivel, 0.4, p, v);
//    }
//};
//
//class CCurbaHilbert
//{
//public:
//    void curbaHilbert(double lungime, int nivel, CPunct& p, CVector& v, int d)
//    {
//        if (nivel == 0)
//        {
//        }
//        else
//        {
//            v.rotatie(d * 90);
//            curbaHilbert(lungime, nivel - 1, p, v, -d);
//
//            v.deseneaza(p, lungime);
//            p = v.getDest(p, lungime);
//
//            v.rotatie(-d * 90);
//            curbaHilbert(lungime, nivel - 1, p, v, d);
//
//            v.deseneaza(p, lungime);
//            p = v.getDest(p, lungime);
//
//            curbaHilbert(lungime, nivel - 1, p, v, d);
//
//            v.rotatie(-d * 90);
//            v.deseneaza(p, lungime);
//            p = v.getDest(p, lungime);
//
//            curbaHilbert(lungime, nivel - 1, p, v, -d);
//
//            v.rotatie(d * 90);
//        }
//    }
//
//    void afisare(double lungime, int nivel)
//    {
//        CVector v(0.0, 1.0);
//        CPunct p(0.0, 0.0);
//
//        curbaHilbert(lungime, nivel, p, v, 1);
//    }
//};
/*
class Square {
public:
    CPunct left_down;
    CPunct left_up;
    CPunct right_up;
    CPunct right_down;
    double lungime;

    Square(CPunct left_down, CPunct left_up, CPunct right_up, CPunct right_down, double lungime)
    {
        this->left_down = left_down;
        this->left_up = left_up;
        this->right_up = right_up;
        this->right_down = right_down;
        this->lungime = lungime;
    }

    Square getMiddleSquare()
    {
        double x, y;
        left_down.getxy(x, y);
        CPunct new_left_down = CPunct(x + lungime / 3, y + lungime / 3);
        CPunct new_left_up = CPunct(x + lungime / 3, y + (lungime / 3) * 2);
        CPunct new_right_up = CPunct(x + (lungime / 3) * 2, y + (lungime / 3) * 2);
        CPunct new_right_down = CPunct(x + (lungime / 3) * 2, y + lungime / 3);

        return Square(new_left_down, new_left_up, new_right_up, new_right_down, lungime / 3);
    }

    Square shiftSquare(Square square, double x_axis, double y_axis)
    {
        double x, y;
        square.left_down.getxy(x, y);
        CPunct new_left_down = CPunct(x + x_axis, y + y_axis);
        square.left_up.getxy(x, y);
        CPunct new_left_up = CPunct(x + x_axis, y + y_axis);
        square.right_up.getxy(x, y);
        CPunct new_right_up = CPunct(x + x_axis, y + y_axis);
        square.right_down.getxy(x, y);
        CPunct new_right_down = CPunct(x + x_axis, y + y_axis);
        return Square(new_left_down, new_left_up, new_right_up, new_right_down, square.lungime);
    }

    vector<Square> getOtherSquares()
    {
        Square middle_square = getMiddleSquare();
        vector<Square> result;

        result.push_back(shiftSquare(middle_square, -middle_square.lungime, middle_square.lungime));
        result.push_back(shiftSquare(middle_square, 0, middle_square.lungime));
        result.push_back(shiftSquare(middle_square, middle_square.lungime, middle_square.lungime));
        result.push_back(shiftSquare(middle_square, middle_square.lungime, 0));
        result.push_back(shiftSquare(middle_square, middle_square.lungime, -middle_square.lungime));
        result.push_back(shiftSquare(middle_square, 0, -middle_square.lungime));
        result.push_back(shiftSquare(middle_square, -middle_square.lungime, -middle_square.lungime));
        result.push_back(shiftSquare(middle_square, -middle_square.lungime, 0));

        return result;
    }

    void draw() 
    {
        glColor3f(1.0, 0.1, 0.1);
        glBegin(GL_LINE_LOOP);
        double x, y;
        left_down.getxy(x, y);
        glVertex2d(x, y);
        left_up.getxy(x, y);
        glVertex2d(x, y);
        right_up.getxy(x, y);
        glVertex2d(x, y);
        right_down.getxy(x, y);
        glVertex2d(x, y);
        glEnd();
        glFlush();
    }

};

class SierpinskiCarpet
{
public:
    void sierpinskiCarpet(Square square, int nivel)
    {
        if (nivel == 0)
        {
            Square middle_square = square.getMiddleSquare();
            middle_square.draw();
        }
        else
        {
            Square middle_square = square.getMiddleSquare();
            middle_square.draw();

            for (Square sub_square : square.getOtherSquares())
            {
                sierpinskiCarpet(sub_square, nivel - 1);
            }
        }
    }

    void draw(int nivel)
    {
        double starting_x = -0.50;
        double starting_y = -0.50;
        double lungime = 1;
        CPunct p1 = CPunct(starting_x, starting_y);
        CPunct p2 = CPunct(starting_x, starting_y + lungime);
        CPunct p3 = CPunct(starting_x + lungime, starting_y + lungime);
        CPunct p4 = CPunct(starting_x + lungime, starting_y);
        Square initial_square = Square(p1, p2, p3, p4, lungime);
        initial_square.draw();
        sierpinskiCarpet(initial_square, nivel);
    }
};

class Segment {
public:
    CPunct starting_point;
    CPunct ending_point;
    double length;
    double angle;

    Segment(CPunct starting_point, double length, double angle)
    {
        this->starting_point = starting_point;
        this->length = length;
        this->angle = 0;
        double x, y;
        starting_point.getxy(x, y);
        ending_point = CPunct(x + length, y);
        rotate(angle);
    }

    void multiply(double mat[2][2],
        double vec[2],
        double res[2])
    {
        int i, j;
        res[0] = 0;
        res[1] = 0;
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
            {
                res[i] += mat[i][j] * vec[j];
            }
        }
    }
    void rotate(double degrees)
    {
        if (degrees == 0)
        {
            return;
        }
        double rotation_matrix[2][2];
        double coordinate_vector[2];
        double x_start, y_start;
        double x_end, y_end;
        starting_point.getxy(x_start, y_start);
        ending_point.getxy(x_end, y_end);
        coordinate_vector[0] = x_end - x_start;
        coordinate_vector[1] = y_end - y_start;
        if (degrees > 0)
        {
            rotation_matrix[0][0] = cos(degrees);
            rotation_matrix[1][1] = cos(degrees);
            rotation_matrix[0][1] = -sin(degrees);
            rotation_matrix[1][0] = sin(degrees);
        }
        else
        {
            rotation_matrix[0][0] = cos(-degrees);
            rotation_matrix[1][1] = cos(-degrees);
            rotation_matrix[0][1] = sin(-degrees);
            rotation_matrix[1][0] = -sin(-degrees);
        }

        double mul_result[2];

        multiply(rotation_matrix, coordinate_vector, mul_result);

        mul_result[0] += x_start;
        mul_result[1] += y_start;

        ending_point = CPunct(mul_result[0], mul_result[1]);
        angle += degrees;
    }

    void draw()
    {
        glColor3f(1.0, 0.1, 0.1);
        glBegin(GL_LINE_STRIP);
        double x, y;
        starting_point.getxy(x, y);
        glVertex2d(x, y);
        ending_point.getxy(x, y);
        glVertex2d(x, y);
        glEnd();
    }
};


class SierpinkiArrowhead {
public:

    CPunct sierpinski_arrowhead(int level, Segment segment, double angle)
    {
        if (level == 0)
        {
            segment.draw();
            return segment.ending_point;
        }
        else
        {
            Segment s1 = Segment(segment.starting_point, segment.length / 2, segment.angle);
            CPunct p1 = sierpinski_arrowhead(level - 1, s1, angle);
            Segment s2 = Segment(p1, segment.length / 2, s1.angle - angle);
            CPunct p2 = sierpinski_arrowhead(level - 1, s2, -angle);
            Segment s3 = Segment(p2, segment.length / 2, s2.angle - angle);
            if (level % 2 == 0)
            {
                s3.rotate(4 * angle);
            }
            return sierpinski_arrowhead(level - 1, s3, angle);

        }
    }
    void display(double length, int level)
    {
        Segment s = Segment(CPunct(-0.5, 0.75), length, -deg_90);

        if (level % 2 == 0)
        {
            sierpinski_arrowhead(level, s, -deg_60);
        }
        else
        {
            s.rotate(deg_60);
            sierpinski_arrowhead(level, s, deg_60);
        }
    }
};


class InvertedPerronTree
{
public:
    void arborePerron(double lungime,
        int nivel,
        double factordiviziune,
        CPunct p,
        CVector v)
    {
        assert(factordiviziune != 0);
        CPunct p1, p2, p3, p4, p5, p6;
        if (nivel == 0)
        {
        }
        else
        {
            v.rotatie(45);
            v.deseneaza(p, lungime);
            p2 = v.getDest(p, lungime);

            v.rotatie(-90);
            v.deseneaza(p, lungime);
            p1 = v.getDest(p, lungime);
            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p1, v);

            v.rotatie(105);
            v.deseneaza(p2, lungime);
            p3 = v.getDest(p2, lungime);
            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p3, v);

            v.rotatie(-60);
            v.deseneaza(p2, lungime);
            p4 = v.getDest(p2, lungime);

            v.rotatie(-90);
            v.deseneaza(p4, lungime * factordiviziune);
            p5 = v.getDest(p4, lungime * factordiviziune);
            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p5, v);

            v.rotatie(120);
            v.deseneaza(p4, lungime * factordiviziune);
            p6 = v.getDest(p4, lungime * factordiviziune);
            arborePerron(lungime * factordiviziune, nivel - 1, factordiviziune, p6, v);
        }
    }

    void afisare(double lungime, int nivel)
    {
        CVector v(0.0, 1.0);
        CPunct p(0.0, 0.8);

        v.rotatie(180);
        v.deseneaza(p, 0.1);
        p = v.getDest(p, 0.1);
        arborePerron(lungime, nivel, 0.4, p, v);
    }
};*/
/*
// afisare curba lui Koch "fulg de zapada"
void Display1() {
    CCurbaKoch cck;
    cck.afisare(sqrt(3.0), nivel);

    char c[3];
    sprintf(c, "%2d", nivel);
    glRasterPos2d(-0.98, -0.98);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);

    glRasterPos2d(-1.0, 0.9);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'c');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'u');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'b');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'a');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'u');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'K');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'o');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'c');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'h');

    nivel++;
}

// afisare arbore binar
void Display2() {
    CArboreBinar cab;
    cab.afisare(1, nivel);

    char c[3];
    sprintf(c, "%2d", nivel);
    glRasterPos2d(-0.98, -0.98);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);

    glRasterPos2d(-1.0, 0.9);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'a');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'b');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'o');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'b');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'n');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'a');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');

    nivel++;
}

// afisare arborele lui Perron
void Display3() {
    CArborePerron cap;

    char c[3];
    sprintf(c, "%2d", nivel);
    glRasterPos2d(-0.98, -0.98);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);

    glRasterPos2d(-1.0, -0.9);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'a');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'b');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'o');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'P');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'o');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'n');

    glPushMatrix();
    glLoadIdentity();
    glScaled(0.4, 0.4, 1);
    glTranslated(-0.5, -0.5, 0.0);
    cap.afisare(1, nivel);
    glPopMatrix();
    nivel++;
}

// afisare curba lui Hilbert
void Display4() {
    CCurbaHilbert cch;
    cch.afisare(0.05, nivel);

    char c[3];
    sprintf(c, "%2d", nivel);
    glRasterPos2d(-0.98, -0.98);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);

    glRasterPos2d(-1.0, -0.9);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'c');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'u');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'b');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'a');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'H');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'b');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'r');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 't');

    nivel++;
}

void Display5()
{
    SierpinskiCarpet sc;
    sc.draw(nivel);
    char c[3];
    sprintf(c, "%2d", nivel);
    glRasterPos2d(-0.98, -0.98);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);

    nivel++;
}

void Display6()
{
    InvertedPerronTree ipt;
    ipt.afisare(0.5, nivel);

    char c[3];
    sprintf(c, "%2d", nivel);
    glRasterPos2d(-0.98, -0.98);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);

    nivel++;
}

void Display7()
{
    SierpinkiArrowhead sa;
    sa.display(1.25, nivel);
    char c[3];
    sprintf(c, "%2d", nivel);
    glRasterPos2d(-0.98, -0.98);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'N');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'i');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'v');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'e');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, 'l');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '=');
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[0]);
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, c[1]);

    nivel++;
}

void Init(void) {

    glClearColor(1.0, 1.0, 1.0, 1.0);

    glLineWidth(1);

    glPointSize(3);

    glPolygonMode(GL_FRONT, GL_LINE);
}

void Display(void)
{
    switch (prevKey)
    {
    case '0':
        glClear(GL_COLOR_BUFFER_BIT);
        nivel = 0;
        fprintf(stderr, "nivel = %d\n", nivel);
        break;
    case '1':
        glClear(GL_COLOR_BUFFER_BIT);
        Display1();
        break;
    case '2':
        glClear(GL_COLOR_BUFFER_BIT);
        Display2();
        break;
    case '3':
        glClear(GL_COLOR_BUFFER_BIT);
        Display3();
        break;
    case '4':
        glClear(GL_COLOR_BUFFER_BIT);
        Display4();
        break;
    case '5':
        glClear(GL_COLOR_BUFFER_BIT);
        Display5();
        break;
    case '6':
        glClear(GL_COLOR_BUFFER_BIT);
        Display6();
        break;
    case '7':
        glClear(GL_COLOR_BUFFER_BIT);
        Display7();
        break;
    case '8':
        glClear(GL_COLOR_BUFFER_BIT);
        Display8();
        break;
    case '9':
        glClear(GL_COLOR_BUFFER_BIT);
        Display9();
        break;
    case 'q':
        glClear(GL_COLOR_BUFFER_BIT);
        DisplayMandelbrot();
        break;
    default:
        break;
    }

    glFlush();
}

void Reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
}

void KeyboardFunc(unsigned char key, int x, int y)
{
    prevKey = key;
    if (key == 27) // escape
        exit(0);
    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {}
*/
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