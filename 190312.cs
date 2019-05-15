using System;
using System.Collections.Generic;
using amateurmathlib;

namespace C_
{
    class Program
    {
        double f1(double x)
        {
            return x * x * Math.Exp(x) - 4;
            
        }

        double f2(double x, double y)
        {
            return x + y;
        }


        static void Main(string[] args)
        {

            
        //выявление ошибок
        /*
        try
        {
            Matrix m1(2, 2);
            m1(0, 0) = 1;
            m1(0, 1) = 2;
            m1(1, 0) = 3;
            m1(1, 1) = 4;

            Matrix m2(2, 2);
            m2(0, 0) = 5;
            m2(0, 1) = 6;
            m2(1, 0) = 7;
            m2(1, 1) = 8;

            Matrix m3 = m1 + m2;

            cout << m3(1, 0) << endl;
        }
        catch (string s)
        {
            cout << s << endl;
        }

        */

           Matrix a = new Matrix(3, 3);

            a[0, 0] = 10;
            a[0, 1] = 2;
            a[0, 2] = 3;

            a[1, 0] = 4;
            a[1, 1] = 10;
            a[1, 2] = 2;

            a[2, 0] = 3;
            a[2, 1] = 4;
            a[2, 2] = 10;



            Matrix b = new Matrix(3, 1);

            b[0] = 2;
            b[1] = 0;
            b[2] = 1;

            List<Matrix.Funcs> funcs = new List<Matrix.Funcs>(2);

            funcs[0] =  (double x, List<double> y) => 
            {
                return y[1];
            };
            funcs[1] = (double x, List<double> y) => 
            {
                return x + 2 * y[0] - y[1];
            };
            double x0 = 0;
            List<double> y0 = new List<double>{ 3.0 / 4.0, 1.0 / 2.0 };
            S_Ode_Solver euler = new S_Ode_Solver(s_ode_euler_koshi, 2);
            List<double> y = euler[funcs, x0, y0, 2.0, 0.0001];

            foreach (var el in y)
            {
                Console.WriteLine(el);
            }

        }
    }
}
