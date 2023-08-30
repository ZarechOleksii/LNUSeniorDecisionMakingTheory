using System;
using System.Windows.Forms;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using System.Collections.Generic;
using System.Windows.Forms.DataVisualization.Charting;
using System.Linq;

namespace Lab5
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        double[,] x1, x2, x3;
        double[,] y;
        double[,] yh;

        int n1, n2, n3, n;
        int size;

        private void CalculateButton_Click(object sender, EventArgs e)
        {
            CleanUp();

            int p1, p2, p3;
            p1 = int.Parse(polPow1.Text);
            p2 = int.Parse(polPow2.Text);
            p3 = int.Parse(polPow3.Text);

            x1 = Utils.Normalization(x1);
            x2 = Utils.Normalization(x2);
            x3 = Utils.Normalization(x3);
            y = Utils.Normalization(y);

            lamdaTextBox.Text += "\nLambdas \n\n";
            double[] L = Lambda(p1, p2, p3);
            for (int i = 0; i < L.Length; ++i)
            {
                lamdaTextBox.Text += L[i].ToString();
                lamdaTextBox.Text += "\n";
            }

            double[,] Lambda1 = new double[n1, p1 + 1];
            double[,] Lambda2 = new double[n2, p2 + 1];
            double[,] Lambda3 = new double[n3, p3 + 1];

            for (int i = 0; i < n1; ++i)
            {
                Lambda1[i, 0] = L[0]; // перший коефіцієнт розкладу для всіх однаковий
                for (int j = 1; j < p1 + 1; ++j)
                {
                    Lambda1[i, j] = L[i * p1 + j];
                }
            }
            for (int i = 0; i < n2; ++i)
            {
                Lambda2[i, 0] = L[0];
                for (int j = 1; j < p2 + 1; ++j)
                {
                    Lambda2[i, j] = L[n1 * p1 + i * p2 + j];
                }
            }
            for (int i = 0; i < n3; ++i)
            {
                Lambda3[i, 0] = L[0];
                for (int j = 1; j < p3 + 1; ++j)
                {
                    Lambda3[i, j] = L[n1 * p1 + n2 * p2 + i * p3 + j];
                }
            }

            double[,] A1 = Utils.GetA(Utils.GetPsi(size, Lambda1, x1), y);
            double[,] A2 = Utils.GetA(Utils.GetPsi(size, Lambda2, x2), y);
            double[,] A3 = Utils.GetA(Utils.GetPsi(size, Lambda3, x3), y);

            matrixATextBox.Text += "\nMatrix A coefficients\n\n";
            for (int i = 0; i < A1.GetLength(0); ++i)
            {
                for (int j = 0; j < A1.GetLength(1); ++j)
                {
                    matrixATextBox.Text += A1[i, j].ToString() + "    ";
                }
                matrixATextBox.Text += "\n";
            }
            matrixATextBox.Text += "\n";

            for (int i = 0; i < A2.GetLength(0); ++i)
            {
                for (int j = 0; j < A2.GetLength(1); ++j)
                {
                    matrixATextBox.Text += A2[i, j].ToString() + "    ";
                }
                matrixATextBox.Text += "\n";
            }
            matrixATextBox.Text += "\n";

            for (int i = 0; i < A3.GetLength(0); ++i)
            {
                for (int j = 0; j < A3.GetLength(1); ++j)
                {
                    matrixATextBox.Text += A3[i, j].ToString() + "    ";
                }
                matrixATextBox.Text += "\n";
            }
            matrixATextBox.Text += "\n";

            matrixCTextBox.Text += "\nMatrix C coefficients\n\n";

            double[][] C1 = new double[4][];
            for (int i = 0; i < n; ++i)
            {
                var psi1 = Utils.GetPsi(size, Lambda1, x1);
                var psi2 = Utils.GetPsi(size, Lambda2, x2);
                var psi3 = Utils.GetPsi(size, Lambda3, x3);
                C1[i] = Utils.GetC(size, y, i, A1, A2, A3, psi1, psi2, psi3);
            }
            for (int i = 0; i < C1.Length; ++i)
            {
                for (int j = 0; j < C1[i].Length; ++j)
                {
                    matrixCTextBox.Text += C1[i][j].ToString() + "    ";
                }
                matrixCTextBox.Text += "\n";
            }

            yh = new double[size, n];
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    double sum = 0;
                    sum += C1[j][0] * Math.Log(Utils.GetFi(size, A1, Utils.GetPsi(size, Lambda1, x1), 0)[i] + 1);
                    sum += C1[j][1] * Math.Log(Utils.GetFi(size, A2, Utils.GetPsi(size, Lambda2, x2), 1)[i] + 1);
                    sum += C1[j][2] * Math.Log(Utils.GetFi(size, A3, Utils.GetPsi(size, Lambda3, x3), 2)[i] + 1);
                    yh[i, j] = Math.Exp(sum) - 1;
                }
            }

            for (int j = 0; j < n; ++j)
            {
                functionsTextBox.Text += "f" + (j + 1) + "(x1, x2, x3)" + " = e^(" + Math.Abs(C1[j][0]).ToString() + " * ln(1 + f" + (j + 1) + "1(x1)" + ")";
                if (C1[j][1] > 0)
                {
                    functionsTextBox.Text += " + \n+";
                }
                else
                {
                    functionsTextBox.Text += " - \n- ";
                }
                functionsTextBox.Text += Math.Abs(C1[j][1]).ToString() + " * ln(1 + f" + (j + 1) + "2(x2)" + ")";
                if (C1[j][2] > 0)
                {
                    functionsTextBox.Text += " + \n+";
                }
                else
                {
                    functionsTextBox.Text += " - \n- ";
                }
                functionsTextBox.Text += Math.Abs(C1[j][2]).ToString() + " * ln(1 + f" + (j + 1) + "3(x3)" + "))\n";
            }

            var charts = new Chart[] { chart1, chart2, chart3, chart4 };
            var errors = new TextBox[] { errorbox1, errorbox2, errorbox3, errorbox4 };
            List<double> values = new List<double>();
            foreach (var series in charts.SelectMany(chart => chart.Series))
            {
                series.Points.Clear();
            }

            for (int m = 0; m < charts.Length; m++)
            {
                double max = 0;
                for (int i = 0; i < size; ++i)
                {
                    
                    charts[m].Series[0].Points.AddXY(i, y[i, m]);
                    charts[m].Series[1].Points.AddXY(i, yh[i, m]);
                    if (Math.Abs(y[i, m] - yh[i, m]) > max)
                    {
                        max = Math.Abs(y[i, m] - yh[i, m]);
                    }
                }
                values.Add(max);
                errors[m].Text = max.ToString();
            }
        }

        private void OpenFile_Click(object sender, EventArgs e)
        {
            openFileDialog1.InitialDirectory = Directory.GetCurrentDirectory();
            if (openFileDialog1.ShowDialog() == DialogResult.Cancel)
                return;

            CleanUp();

            string fileName = openFileDialog1.FileName;

            // Initialize inputs
            n1 = int.Parse(textBox1.Text);
            n2 = int.Parse(textBox2.Text);
            n3 = int.Parse(textBox3.Text);
            n = int.Parse(textBox4.Text);
            size = int.Parse(textBox8.Text);

            x1 = new double[size, n1];
            x2 = new double[size, n2];
            x3 = new double[size, n3];
            y = new double[size, n];

            //Read from file
            string[] s = File.ReadAllLines(fileName);
            if (s.Length != x1.GetLength(0) + 1)
            {
                throw new Exception("Dimensions do not much !");
            }

            for (int i = 0; i < s.Length - 1; ++i)
            {
                string[] arr_s = s[i + 1].Split(';');
                x1[i, 0] = Convert.ToDouble(arr_s[0]);
                x1[i, 1] = Convert.ToDouble(arr_s[1]);
                x2[i, 0] = Convert.ToDouble(arr_s[2]);
                x2[i, 1] = Convert.ToDouble(arr_s[3]);
                x3[i, 0] = Convert.ToDouble(arr_s[4]);
                x3[i, 1] = Convert.ToDouble(arr_s[5]);
                x3[i, 2] = Convert.ToDouble(arr_s[6]);
                y[i, 0] = Convert.ToDouble(arr_s[7]);
                y[i, 1] = Convert.ToDouble(arr_s[8]);
                y[i, 2] = Convert.ToDouble(arr_s[9]);
                y[i, 3] = Convert.ToDouble(arr_s[10]);
            }
        }

        private double[] Lambda(int p1, int p2, int p3)
        {
            double[] matrix_b = Utils.GetB(y);

            int m = n1 * p1 + n2 * p2 + n3 * p3 + 1;
            var matrix = Matrix<double>.Build.Dense(size, m);

            for (int i = 0; i < size; ++i)
            {
                matrix[i, 0] = Math.Log(1 + 0.5); // 1 + T0
                for (int j = 0; j < n1; ++j)
                {
                    for (int p = 1; p < p1 + 1; ++p)
                    {
                        matrix[i, j * p1 + p] = Math.Log(1e-9 + 1 + Utils.GetTStar(x1[i, j], p));
                    }
                }
                for (int j = 0; j < n2; ++j)
                {
                    for (int p = 1; p < p2 + 1; ++p)
                    {
                        matrix[i, n1 * p1 + j * p2 + p] = Math.Log(1e-9 + 1 + Utils.GetTStar(x2[i, j], p));
                    }
                }
                for (int j = 0; j < n3; ++j)
                {
                    for (int p = 1; p < p3 + 1; ++p)
                    {
                        matrix[i, n1 * p1 + n2 * p2 + j * p3 + p] = Math.Log(1e-9 + 1 + Utils.GetTStar(x3[i, j], p));
                    }
                }
            }

            var vector = Matrix<double>.Build.Dense(size, 1);
            for (int i = 0; i < size; ++i)
            {
                vector[i, 0] = Math.Log(1e-9 + 1 + matrix_b[i]);
            }

            var resVector = Utils.SmallerSquareMethod(matrix, vector);

            double[] res = new double[m];
            for (int i = 0; i < m; ++i)
            {
                res[i] = resVector[i, 0];
            }

            return res;
        }

        private void CleanUp()
        {
            lamdaTextBox.Text = "";
            matrixATextBox.Text = "";
            matrixCTextBox.Text = "";
            functionsTextBox.Text = "";
        }
    }
}
