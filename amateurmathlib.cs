using System;
using System.Collections.Generic;

namespace amateurmathlib
{
	class Matrix
	{
		//передаем функцию
	 	public delegate double Func(double a);
		public delegate double Func2(double a, double b);
		public delegate double Funcs(double a, List<double> list);

		public string ToString(Matrix m) //вывод матрицы
		{
			string str = string.Empty;
			for (int i = 0; i < m.nRows(); i++)
			{
				for (int j = 0; j < m.nColumns(); j++)
				{
					double tmp;										//ммм, костыли
					tmp = m[i, j];
					if (Math.Abs(tmp) < 0.000001)
						tmp = 0.0;
					str += (Console.WindowWidth = 9)  + "" + tmp + " \t";
				}
				str += "\n";
			}
			return str;
		}


		public double  this[int row, int column = 0]					//выражаем матрицу в виде вектора
		{
			/*
			234
			567
			108

			(2,3,4,5,6,7,1,0,8)
			*/
			get
			{
				if (row < 0 || row >= n_rows)
					throw new Exception("Wrong row index");
				if (column < 0 || column >= n_columns)
					throw new Exception("Wrong column index");
				return el[n_columns * row + column];
			}
			set
			{
				el[row] = value;
			}
		}

		public Matrix(int n_rows, int n_columns = 1)											//объявление
		{
			if (n_rows < 1 || n_columns < 1)
				throw new Exception("Wrong number of rows and columns");						//выкинуть ошибку
																					//this-> указатель на объект класса
			this.n_rows = n_rows;
			this.n_columns = n_columns;
		 	el = new List<double>(n_rows * n_columns);
		}
		private	int nRows()
		{
			return n_rows;
		}
		private	int nColumns()
		{
			return n_columns;
		}

		private static Matrix add(Matrix m1, Matrix m2, int sign)						//костыль для сложения-вычитания
		{
			Matrix result = new Matrix(m1.nRows(), m1.nColumns());
			for (int i = 0; i < m1.nRows(); i++)
			{
				for (int j = 0; j < m1.nColumns(); j++)
				{
					result[i, j] = m1[i, j] + sign * m2[i, j];
				}
			}
			return result;
		}

		public static  Matrix operator +(Matrix m1, Matrix m2)
		{
			return add(m1, m2, 1);
		}

		public static  Matrix operator -(Matrix m1, Matrix m2)
		{
			return add(m1, m2, -1);
		}

		public static Matrix operator *(Matrix m, double value)						//умножение матрицы на константу
		{
			Matrix result = new Matrix(m.nRows(), m.nColumns());
			for (int i = 0; i < m.nRows(); i++)
			{
				for (int j = 0; j < m.nColumns(); j++)
				{
					result[i, j] = value * m[i, j];
				}
			}
			return result;
		}

		public static Matrix operator *(double value, Matrix m)						//умножение константы на матрицу, для идиотов
		{
			return m * value;
		}

		public static Matrix operator *(Matrix m1, Matrix m2)							//умножение матрицы на матрицу
		{
			Matrix result = new Matrix(m1.nRows(), m2.nColumns());
			for (int i = 0; i < m1.nRows(); i++)
			{

				for (int j = 0; j < m2.nColumns(); j++)
				{
					double tmp = 0;
					for (int k = 0; k < m2.nRows(); k++)
					{
						tmp += m1[i, k] * m2[k, j];
					}
					result[i, j] = tmp;
				}
			}
			return result;
		}



		double det(Matrix m)
		{
			for (int i = 0; i < m.nRows() - 1; i++)
			{
				for (int j = i + 1; j < m.nRows(); j++)
				{
					double coeff = m[j, i] / m[i, i];
					for (int k = i; k < m.nColumns(); k++)
					{
						m[j, k] = m[j, k] - coeff * m[i, k];
					}
				}
			}

			double result = 1;
			for (int i = 0; i < m.nRows(); i++)
			{
				result *= m[i, i];
			}
			return result;
		}

		Matrix jordan_gauss(Matrix a, Matrix b)
		{
			for (int i = 0; i < a.nRows() - 1; i++)					//ведущая строка
			{
				for (int j = i + 1; j < a.nRows(); j++)				//текущая строка
				{
					double coeff = a[j, i] / a[i, i];				//коэффициент
					for (int k = i; k < a.nColumns(); k++)
					{
						a[j, k] = a[j, k] - coeff * a[i, k];
					}
					b[j, 0] = b[j, 0] - coeff * b[i, 0];
				}


			}

			for (int i = a.nRows() - 1; i >= 1; i--)
			{
				for (int j = i - 1; j >= 0; j--)
				{
					double coeff = a[j, i] / a[i, i];
					a[j, i] = a[j, i] - coeff * a[i, i];
					b[j, 0] = b[j, 0] - coeff * b[i, 0];
				}
			}

			for (int i = 0; i < a.nRows(); i++)
			{
				b[i, 0] = b[i, 0] / a[i, i];
				a[i, i] = 1;
			}

			return b;
		}

		Matrix inv_jg(Matrix a)
		{

			if (a.nRows() != a.nColumns())
				throw new Exception("Matrix is not square");

			if (Math.Abs(det(a)) < 0.00000000000001)										//вещественное число
				throw new Exception(" Det cannot be equal to 0");

			Matrix result = new Matrix(a.nRows(), a.nColumns());
			for (int i = 0; i < a.nRows(); i++)
			{
				result[i, i] = 1;													//единичная матрица
			}

			for (int i = 0; i < a.nRows() - 1; i++)					//ведущая строка
			{
				for (int j = i + 1; j < a.nRows(); j++)				//текущая строка
				{
					double coeff = a[j, i] / a[i, i];				//коэффициент
					for (int k = i; k < a.nColumns(); k++)
					{
						a[j, k] = a[j, k] - coeff * a[i, k];
					}
					for (int k = 0; k < a.nColumns(); k++)
					{
						result[j, k] = result[j, k] - coeff * result[i, k];
					}
				}
			}

			for (int i = a.nRows() - 1; i >= 1; i--)
			{
				for (int j = i - 1; j >= 0; j--)
				{
					double coeff = a[j, i] / a[i, i];
					a[j, i] = a[j, i] - coeff * a[i, i];
					for (int k = 0; k < a.nColumns(); k++)
					{
						result[j, k] = result[j, k] - coeff * result[i, k];
					}

				}
			}

			for (int i = 0; i < a.nRows(); i++)
			{
				for (int k = 0; k < a.nColumns(); k++)
				{
					result[i, k] /= a[i, i];
				}
				a[i, i] = 1;
			}

			return result;

		}

		double vect_norm1(Matrix v)
		{
			double result = 0;
			for (int i = 0; i < v.nRows(); i++)
			{
				if ( Math.Abs( v[i] ) > result )
				{
					result = Math.Abs( v[i] );
				}
			}

			return result;
		}

		double vect_norm2(Matrix v)
		{
			double result = 0;
			for (int i = 0; i < v.nRows(); i++)
			{
				result += Math.Abs( v[i] );
			}
			return result;
		}

		double vect_norm3(Matrix v)
		{
			double result = 0;
			for (int i = 0; i < v.nRows(); i++)
			{
				result += Math.Pow(v[i], 2);
			}
			return Math.Sqrt(result);
		}

		double matr_norm1(Matrix m)
		{
			double result = 0;
			for (int i = 0; i < m.nRows(); i++)
			{
				double sum = 0;
				for (int j = 0; j < m.nColumns(); j++)
				{
					sum += Math.Abs( m[i,j] );
				}
				if (result < sum)
				{
					result = sum;
				}
			}
			return result;
		}

		double matr_norm2(Matrix m)
		{
			double result = 0;
			for (int j = 0; j < m.nColumns(); j++)
			{
				double sum = 0;
				for (int i = 0; i < m.nRows(); i++)
				{
					sum += Math.Abs(m[i, j]);
				}
				if (result < sum)
				{
					result = sum;
				}
			}
			return result;
		}


		Matrix jacobi_iter(Matrix m, Matrix b, double eps)
		{
			Matrix result = b;
			for (int k = 0; k < m.nColumns(); k++)
			{
				b[k] /= m[k, k];
			}

			for (int i = 0; i < m.nRows(); i++)
			{
				for (int k = 0; k < m.nColumns(); k++)
				{
					if (k == i)
					{
						continue;
					}

					m[i, k] /= -m[i, i];
				}

				m[i, i] = 0;

			}

			while (true)
			{
				Matrix last = result;
				result = m * last + b;							//метод простых итерация
				if (vect_norm1(result - last) < eps)
				{
					break;
				}
			}


			return result;
		}

		double newton(Func f, double x, double eps)
		{
			while (true)
			{
				double last = x;
				double df = (f(x + eps) - f(x)) / eps;
				x = last - f(x) / df;
				if ( Math.Abs(x-last) < eps )
				{
					break;
				}
			}
			return x;
		}

		double m_hord(Func f, double x, double eps)
		{
			double a = x - eps; 
			while (true)
			{
				double last = x;
				x = last - f(last) * (a - last) / ( f(a) - f(last) );
				if (Math.Abs(x - last) < eps)
				{
					break;
				}
			}
			return x;
		}

		double m_iter(Func f, double x, double eps)
		{
			double a = x - eps;
			while (true)
			{
				double last = x;

				x = last + f(x);

				if (Math.Abs(x - last) < eps)
				{
					break;
				}
			}
			return x;
		}

		double helper_int_rect(Func f, double a, double b, int n)
		{
			double h = (b - a) / n;
			double sum = 0;
			for (double x = a; x < b; x += h)
			{
				sum += f(x) * h;
			}

			return sum;

		}

		double int_rect(Func f, double a, double b, double eps)		//интегрирование методом прямоугольников
		{
			int n = 2;
			double I = helper_int_rect(f, a, b, n);
			while (true)
			{
				
				n *= 2;
				double I_half = helper_int_rect(f, a, b, n);
				double R = (I - I_half) / 1;
				if (Math.Abs(R) < eps)
				{
					return I_half;
				}
				I = I_half;
			}
		}

		double helper_int_trap(Func f, double a, double b, int n)
		{
			double h = (b - a) / n;
			double sum = ( f(a) + f(b) ) / 2.0;
			for (double x = a+h; x < b-h; x += h)
			{
				sum += f(x);
			}

			return sum * h;

		}

		double int_trap(Func f, double a, double b, double eps)		
		{
			int n = 2;
			double I = helper_int_trap(f, a, b, n);
			while (true)
			{

				n *= 2;
				double I_half = helper_int_trap(f, a, b, n);
				double R = (I_half - I) / 3.0;
				if (Math.Abs(R) < eps)
				{
					return I_half + R;
				}
				I = I_half;
			}
		}

		double helper_int_simpson(Func f, double a, double b, double n)
		{
			double h = (b - a) / n;		//на какие интервалы делится отрезок
			double double_h = 2 * h;

			double result_1 = f(a) + f(b);

			double result_2 = 0.0;
			for (double x = a + h; x < b; x += double_h)
			{
				result_2 += f(x);
			}

			double result_3 = 0.0;
			for (double x = a + double_h; x < b; x += double_h)
			{
				result_3 += f(x);
			}

			return h * (result_1 + 4 * result_2 + 2 * result_3) / 3.0;
		}

		double int_simpson(Func f, double a, double b, double eps)
		{
			int n = 2;
			double I = helper_int_simpson(f, a, b, n);
			while (true)
			{

				n *= 2;
				double I_half = helper_int_simpson(f, a, b, n);
				double R = (I_half - I) / 15.0;
				if (Math.Abs(R) < eps)
				{
					return I_half + R;
				}
				I = I_half;
			}
		}

		List<double> s_ode_euler_koshi(List<Funcs> f, double x, List<double> y, double h)
		{
			List<double> result = new List<double>(y.Count);
			for (int k = 0; k < y.Count; k++)
			{
				result[k] = y[k] + h * f[k](x, y);
			}
			double h_half = h / 2;
			for (int k = 0; k < y.Count; k++)
			{
				result[k] = y[k] + h_half * ( f[k](x,y) + f[k](x+h, result));
			}
			return result;
		} 


		private	List<double> el;
		private	int n_rows;																	//строки
		private	int n_columns;																//столбцы
	};

	class S_Ode_Solver
	{
		delegate double Funcs(double a, List<double> list);

		delegate List<double> Method(List<Funcs> lst1, double a, List<double> lst, double b);
		public S_Ode_Solver(Method method, int order) //конструктор
		{
			this.method = method;
			this.order = order;
		}

		public List<double> this[List<Matrix.Funcs> f, double x0, List<double> y0, double x, double eps]//выражаем матрицу в виде вектора
		{
			get
			{
				int n = 2;
				List<double> y = solve(f, x0, y0, x, n);
				double magic_denom = (Math.Pow(2, order) - 1);
				while (true)
				{

					n *= 2;
					List<double> y_half = solve(f, x0, y0, x, n);

					double max_R = 0;
					for (int k = 0; k < y.Count; k++)
					{
						double R = Math.Abs(y_half[k] - y[k])/magic_denom;
						if (R > max_R)
						{
							max_R = R;
						}
					}
					//cout << max_R << endl;
					if (max_R < eps)
					{
						for (int i = 0; i < y.Count; i++)
						{
							y_half[i] += ((y_half[i] - y[i]) / (magic_denom));
						}
						return y_half;
					}
					y = y_half;
				}
			}
		}

		private List<double> solve(List<Funcs> f, double x0, List<double> y0, double x, int n)
		{
			List<double> y_i = y0;
			double h = (x - x0) / n;


			for (double x_i = x0; x_i <= x; x_i += h)
			{
				y_i = method(f, x_i, y_i, h);
			}

			return y_i;
		}
		private	Method method;
		private	int order;
	}
};