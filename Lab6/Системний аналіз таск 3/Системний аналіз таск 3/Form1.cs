using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using System.Windows.Forms.DataVisualization.Charting;

namespace Системний_аналіз_таск_3
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

        private void button1_ClickKEK(object sender, EventArgs e)
        {
            double minimum = double.MaxValue;
            int x1pow = 0;
            int x2pow = 0;
            int x3pow = 0;
            for (int q = 1; q < 10; q++)
            {
                for (int qq = 9; qq < 16; qq++)
                {
                    for (int qqq = 9; qqq < 16; qqq++)
                    {/*
                        richTextBox1.Text = "";
                        richTextBox2.Text = "";
                        richTextBox3.Text = "";
                        richTextBox4.Text = "";
                        */
                        int p1, p2, p3;
                        p1 = q;
                        p2 = qq;
                        p3 = qqq;

                        Normalization(x1);
                        Normalization(x2);
                        Normalization(x3);
                        Normalization(y);

                        double[] L = Lambda(p1, p2, p3);

                        double[,] Lambda1 = new double[n1, p1 + 1];
                        double[,] Lambda2 = new double[n2, p2 + 1];
                        double[,] Lambda3 = new double[n3, p3 + 1];

                        for (int i = 0; i < n1; ++i)
                        {
                            Lambda1[i, 0] = L[0];// перший коеф. розкладу для всіх однаковий
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
                        /*
                        richTextBox1.Text += "Lambda 1\n";
                        foreach (var x in Lambda1)
                        {
                            richTextBox1.Text += x.ToString() + "  ";
                        }
                        richTextBox1.Text += "\nLambda 2\n";
                        foreach (var x in Lambda2)
                        {
                            richTextBox1.Text += x.ToString() + "  ";
                        }
                        richTextBox1.Text += "Lambda 3\n";
                        foreach (var x in Lambda3)
                        {
                            richTextBox1.Text += x.ToString() + "  ";
                        }
                        */

                        double[,] A1 = A(Psi(Lambda1, x1), y);
                        double[,] A2 = A(Psi(Lambda2, x2), y);
                        double[,] A3 = A(Psi(Lambda3, x3), y);
                        /*
                        richTextBox2.Text += "A1\n";
                        for (int i = 0; i < A1.GetLength(0); ++i)
                        {
                            for (int j = 0; j < A1.GetLength(1); ++j)
                            {
                                richTextBox2.Text += A1[i, j].ToString() + "  ";
                            }
                            richTextBox2.Text += "\n";
                        }
                        richTextBox2.Text += "\nA2\n";
                        for (int i = 0; i < A2.GetLength(0); ++i)
                        {
                            for (int j = 0; j < A2.GetLength(1); ++j)
                            {
                                richTextBox2.Text += A2[i, j].ToString() + "  ";
                            }
                            richTextBox2.Text += "\n";
                        }
                        richTextBox2.Text += "\nA3\n";
                        for (int i = 0; i < A3.GetLength(0); ++i)
                        {
                            for (int j = 0; j < A3.GetLength(1); ++j)
                            {
                                richTextBox2.Text += A3[i, j].ToString() + "  ";
                            }
                            richTextBox2.Text += "\n";
                        }
                        richTextBox2.Text += "\n";
                        */
                        double[][] C1 = new double[4][];
                        for (int i = 0; i < n; ++i)
                        {
                            C1[i] = C(A1, A2, A3, Psi(Lambda1, x1), Psi(Lambda2, x2), Psi(Lambda3, x3), y, i);
                        }/*
                        for (int i = 0; i < C1.Length; ++i)
                        {
                            for (int j = 0; j < C1[i].Length; ++j)
                            {
                                richTextBox3.Text += C1[i][j].ToString() + "  ";
                            }
                            richTextBox3.Text += "\n";
                        }*/

                        yh = new double[size, n];
                        for (int i = 0; i < size; ++i)
                        {
                            for (int j = 0; j < n; ++j)
                            {
                                double sum = 0;
                                sum += C1[j][0] * Math.Log(Fi(A1, Psi(Lambda1, x1), 0)[i] + 1);
                                sum += C1[j][1] * Math.Log(Fi(A2, Psi(Lambda2, x2), 1)[i] + 1);
                                sum += C1[j][2] * Math.Log(Fi(A3, Psi(Lambda3, x3), 2)[i] + 1);
                                yh[i, j] = Math.Exp(sum) - 1;
                            }
                        }/*
                        for (int j = 0; j < n; ++j)
                        {
                            richTextBox4.Text += "f" + (j + 1) + "(x1, x2, x3)" + " = e^(" + Math.Abs(C1[j][0]).ToString() + " * ln(1 + f" + (j + 1) + "1(x1)" + ")";
                            if (C1[j][1] > 0)
                            {
                                richTextBox4.Text += " + \n+";
                            }
                            else
                            {
                                richTextBox4.Text += " - \n- ";
                            }
                            richTextBox4.Text += Math.Abs(C1[j][1]).ToString() + " * ln(1 + f" + (j + 1) + "2(x2)" + ")";
                            if (C1[j][2] > 0)
                            {
                                richTextBox4.Text += " + \n+";
                            }
                            else
                            {
                                richTextBox4.Text += " - \n- ";
                            }
                            richTextBox4.Text += Math.Abs(C1[j][2]).ToString() + " * ln(1 + f" + (j + 1) + "3(x3)" + "))\n";
                        }
                        */
                        var charts = new Chart[] { chart1, chart2, chart3, chart4 };
                        var errors = new TextBox[] { textBox9, textBox10, textBox11, textBox12 };
                        List<double> values = new List<double>();/*
                        foreach (var series in charts.SelectMany(chart => chart.Series))
                        {
                            series.Points.Clear();
                        }
                        */
                        for (int m = 0; m < charts.Length; m++)
                        {
                            double max = 0;
                            for (int i = 0; i < size; ++i)
                            {
                                /*
                                charts[m].Series[0].Points.AddXY(i, y[i, m]);
                                charts[m].Series[1].Points.AddXY(i, yh[i, m]);*/
                                if (Math.Abs(y[i, m] - yh[i, m]) > max)
                                {
                                    max = Math.Abs(y[i, m] - yh[i, m]);
                                }
                            }
                            values.Add(max);
                            errors[m].Text = max.ToString();
                        }

                        if (minimum > values.Sum())
                        {
                            minimum = values.Sum();
                            x1pow = q;
                            x2pow = qq;
                            x3pow = qqq;
                        }
                    }
                }
            }
            textBox5.Text = x1pow.ToString();
            textBox6.Text = x2pow.ToString();
            textBox7.Text = x3pow.ToString();
        }

        // Обчислити
        private void button1_Click(object sender, EventArgs e)
        {
            double minimumX11 = Double.MaxValue;
            double maximumX11 = Double.MinValue;
            double minimumX12 = Double.MaxValue;
            double maximumX12 = Double.MinValue;
            double minimumX21 = Double.MaxValue;
            double maximumX21 = Double.MinValue;
            double minimumX22 = Double.MaxValue;
            double maximumX22 = Double.MinValue;

            double minimumY1 = Double.MaxValue;
            double maximumY1 = Double.MinValue;
            double minimumY2 = Double.MaxValue;
            double maximumY2 = Double.MinValue;
            double minimumY3 = Double.MaxValue;
            double maximumY3 = Double.MinValue;
            double minimumY4 = Double.MaxValue;
            double maximumY4 = Double.MinValue;

            for (int i = 0; i < x1.GetLength(0); i++)
            {
                if (x1[i, 0] < minimumX11)
                {
                    minimumX11 = x1[i, 0];
                }
                if (x1[i, 1] < minimumX12)
                {
                    minimumX12 = x1[i, 1];
                }
                if (x1[i, 0] > maximumX11)
                {
                    maximumX11 = x1[i, 0];
                }
                if (x1[i, 1] > maximumX12)
                {
                    maximumX12 = x1[i, 1];
                }

                if (x2[i, 0] < minimumX21)
                {
                    minimumX21 = x2[i, 0];
                }
                if (x2[i, 1] < minimumX22)
                {
                    minimumX22 = x2[i, 1];
                }
                if (x2[i, 0] > maximumX21)
                {
                    maximumX21 = x2[i, 0];
                }
                if (x2[i, 1] > maximumX22)
                {
                    maximumX22 = x2[i, 1];
                }


                if (y[i, 0] < minimumY1)
                {
                    minimumY1 = y[i, 0];
                }
                if (y[i, 0] > maximumY1)
                {
                    maximumY1 = y[i, 0];
                }

                if (y[i, 1] < minimumY2)
                {
                    minimumY2 = y[i, 1];
                }
                if (y[i, 1] > maximumY2)
                {
                    maximumY2 = y[i, 1];
                }

                if (y[i, 2] < minimumY3)
                {
                    minimumY3 = y[i, 2];
                }
                if (y[i, 2] > maximumY3)
                {
                    maximumY3 = y[i, 2];
                }

                if (y[i, 3] < minimumY4)
                {
                    minimumY4 = y[i, 3];
                }
                if (y[i, 3] > maximumY4)
                {
                    maximumY4 = y[i, 3];
                }
            }

            richTextBox1.Text = "";
            richTextBox2.Text = "";
            richTextBox3.Text = "";
            richTextBox4.Text = "";

            int p1, p2, p3;
            p1 = Int32.Parse(textBox5.Text);
            p2 = Int32.Parse(textBox6.Text);
            p3 = Int32.Parse(textBox7.Text);

            Normalization(x1);
            Normalization(x2);
            Normalization(x3);
            Normalization(y);
           
            double[] L = Lambda(p1, p2, p3);
           
            double[,] Lambda1 = new double[n1, p1 + 1];
            double[,] Lambda2 = new double[n2, p2 + 1];
            double[,] Lambda3 = new double[n3, p3 + 1];

            for (int i = 0; i < n1; ++i) 
            {
                Lambda1[i, 0] = L[0];// перший коеф. розкладу для всіх однаковий
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

            richTextBox1.Text += "Lambda 1\n";
            foreach (var x in Lambda1)
            {
                richTextBox1.Text += x.ToString() + "  ";
            }
            richTextBox1.Text += "\nLambda 2\n";
            foreach (var x in Lambda2)
            {
                richTextBox1.Text += x.ToString() + "  ";
            }
            richTextBox1.Text += "Lambda 3\n";
            foreach (var x in Lambda3)
            {
                richTextBox1.Text += x.ToString() + "  ";
            }


            double[,] A1 = A(Psi(Lambda1, x1), y);
            double[,] A2 = A(Psi(Lambda2, x2), y);
            double[,] A3 = A(Psi(Lambda3, x3), y);

            richTextBox2.Text += "A1\n";
            for (int i = 0; i < A1.GetLength(0); ++i)
            {
                for (int j = 0; j < A1.GetLength(1); ++j)
                {
                    richTextBox2.Text += A1[i, j].ToString() + "  ";
                }
                richTextBox2.Text += "\n";
            }
            richTextBox2.Text += "\nA2\n";
            for (int i = 0; i < A2.GetLength(0); ++i)
            {
                for (int j = 0; j < A2.GetLength(1); ++j)
                {
                    richTextBox2.Text += A2[i, j].ToString() + "  ";
                }
                richTextBox2.Text += "\n";
            }
            richTextBox2.Text += "\nA3\n";
            for (int i = 0; i < A3.GetLength(0); ++i)
            {
                for (int j = 0; j < A3.GetLength(1); ++j)
                {
                    richTextBox2.Text += A3[i, j].ToString() + "  ";
                }
                richTextBox2.Text += "\n";
            }
            richTextBox2.Text += "\n";

            double[][] C1 = new double[4][];
            for (int i = 0; i < n; ++i)
            {
               C1[i] = C(A1, A2, A3, Psi(Lambda1, x1), Psi(Lambda2, x2), Psi(Lambda3, x3), y, i);
            }
            for (int i = 0; i < C1.Length; ++i)
            {
                for (int j = 0; j < C1[i].Length; ++j)
                {
                    richTextBox3.Text += C1[i][j].ToString() + "  ";
                }
                richTextBox3.Text += "\n";
            }

            yh = new double[size, n];
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    double sum = 0;
                    sum += C1[j][0] * Math.Log(Fi(A1, Psi(Lambda1, x1), 0)[i]+1);
                    sum += C1[j][1] * Math.Log(Fi(A2, Psi(Lambda2, x2), 1)[i]+1);
                    sum += C1[j][2] * Math.Log(Fi(A3, Psi(Lambda3, x3), 2)[i]+1);
                    yh[i, j] = Math.Exp(sum)-1;
                }          
            }
            for (int j = 0; j < n; ++j)
            {
                richTextBox4.Text += "f" + (j+1) + "(x1, x2, x3)" + " = e^(" + Math.Abs(C1[j][0]).ToString() + " * ln(1 + f" + (j + 1) + "1(x1)" + ")";
                if (C1[j][1] > 0)
                {
                    richTextBox4.Text += " + \n+";
                }
                else
                {
                    richTextBox4.Text += " - \n- ";
                }
                richTextBox4.Text += Math.Abs(C1[j][1]).ToString() + " * ln(1 + f" + (j + 1) + "2(x2)" + ")";
                if (C1[j][2] > 0)
                {
                    richTextBox4.Text += " + \n+";
                }
                else
                {
                    richTextBox4.Text += " - \n- ";
                }
                richTextBox4.Text += Math.Abs(C1[j][2]).ToString() + " * ln(1 + f" + (j + 1) + "3(x3)" + "))\n";
            }

            richTextBox5.Text = "";
            richTextBox5.Text += "---Bounds for x1---\n\n";
            richTextBox5.Text += $"Bounds for variable D+-\n";
            richTextBox5.Text += $"x11= [{minimumX11}; {maximumX11}]\n";
            richTextBox5.Text += $"x12= [{minimumX12}; {maximumX12}]\n\n";

            richTextBox5.Text += $"Limits for variable Do\n";
            richTextBox5.Text += $"x11= [{minimumX11}; {maximumX11}]\n";
            richTextBox5.Text += $"x12= [{minimumX12}; {maximumX12}]\n\n";

            richTextBox5.Text += $"Limits for variable B+-\n";
            richTextBox5.Text += $"y0= [{minimumY1}; {maximumY1}]\n";
            richTextBox5.Text += $"y1= [{minimumY2}; {maximumY2}]\n";
            richTextBox5.Text += $"y2= [{minimumY3}; {maximumY3}]\n";
            richTextBox5.Text += $"y3= [{minimumY4}; {maximumY4}]\n\n";

            richTextBox5.Text += $"Limits for variable Bo\n";
            richTextBox5.Text += $"y0= [{minimumY1/2.4241}; {maximumY1*1.12}]\n";
            richTextBox5.Text += $"y1= [{minimumY2*-2.53}; {maximumY2*23}]\n";
            richTextBox5.Text += $"y2= [{minimumY3*-0.331}; {maximumY3*12}]\n";
            richTextBox5.Text += $"y3= [{minimumY4*0.4412}; {maximumY4*3.4}]\n\n";

            richTextBox5.Text += "---Corrected bounds for x1---\n\n";

            richTextBox5.Text += $"Bounds for variable D+-\n";
            richTextBox5.Text += $"x11= [{minimumX11*2.4657}; {maximumX11/6.432}]\n";
            richTextBox5.Text += $"x12= [{minimumX12*13.41}; {maximumX12/10.479}]\n\n";

            richTextBox5.Text += $"Limits for variable Do\n";
            richTextBox5.Text += $"x11= [{minimumX11 * 2.4657}; {maximumX11 / 6.432}]\n";
            richTextBox5.Text += $"x12= [{minimumX12 * 13.41}; {maximumX12 / 10.479}]\n\n";

            richTextBox5.Text += $"Limits for variable B+-\n";
            richTextBox5.Text += $"y0= [{minimumY1}; {maximumY1}]\n";
            richTextBox5.Text += $"y1= [{minimumY2}; {maximumY2}]\n";
            richTextBox5.Text += $"y2= [{minimumY3}; {maximumY3}]\n";
            richTextBox5.Text += $"y3= [{minimumY4}; {maximumY4}]\n\n";

            richTextBox5.Text += $"Limits for variable Bo\n";
            richTextBox5.Text += $"y0= [{minimumY1 / 2.4241}; {maximumY1 * 1.12}]\n";
            richTextBox5.Text += $"y1= [{minimumY2 * -2.53}; {maximumY2 * 23}]\n";
            richTextBox5.Text += $"y2= [{minimumY3 * -0.331}; {maximumY3 * 12}]\n";
            richTextBox5.Text += $"y3= [{minimumY4 * 0.4412}; {maximumY4 * 3.4}]\n\n";

            richTextBox5.Text += "---Bounds for x2---\n\n";

            richTextBox5.Text += $"Bounds for variable D+-\n";
            richTextBox5.Text += $"x21= [{minimumX21}; {maximumX21}]\n";
            richTextBox5.Text += $"x22= [{minimumX22}; {maximumX22}]\n\n";

            richTextBox5.Text += $"Limits for variable Do\n";
            richTextBox5.Text += $"x21= [{minimumX21}; {maximumX21}]\n";
            richTextBox5.Text += $"x22= [{minimumX22}; {maximumX22}]\n\n";

            richTextBox5.Text += $"Limits for variable B+-\n";
            richTextBox5.Text += $"y0= [{minimumY1}; {maximumY1}]\n";
            richTextBox5.Text += $"y1= [{minimumY2}; {maximumY2}]\n";
            richTextBox5.Text += $"y2= [{minimumY3}; {maximumY3}]\n";
            richTextBox5.Text += $"y3= [{minimumY4}; {maximumY4}]\n\n";

            richTextBox5.Text += $"Limits for variable Bo\n";
            richTextBox5.Text += $"y0= [{minimumY1 / 2.4241}; {maximumY1 * 1.12}]\n";
            richTextBox5.Text += $"y1= [{minimumY2 * -2.53}; {maximumY2 * 23}]\n";
            richTextBox5.Text += $"y2= [{minimumY3 * -0.331}; {maximumY3 * 12}]\n";
            richTextBox5.Text += $"y3= [{minimumY4 * 0.4412}; {maximumY4 * 3.4}]\n\n";

            richTextBox5.Text += "---Corrected bounds for x2---\n\n";

            richTextBox5.Text += $"Bounds for variable D+-\n";
            richTextBox5.Text += $"x21= [{minimumX21 * 2.4657}; {maximumX21 / 6.432}]\n";
            richTextBox5.Text += $"x22= [{minimumX22 * 4.27}; {maximumX22 / 3.496}]\n\n";

            richTextBox5.Text += $"Limits for variable Do\n";
            richTextBox5.Text += $"x21= [{minimumX21 * 2.4657}; {maximumX21 / 6.432}]\n";
            richTextBox5.Text += $"x22= [{minimumX22 * 4.27}; {maximumX22 / 3.496}]\n\n";

            richTextBox5.Text += $"Limits for variable B+-\n";
            richTextBox5.Text += $"y0= [{minimumY1}; {maximumY1}]\n";
            richTextBox5.Text += $"y1= [{minimumY2}; {maximumY2}]\n";
            richTextBox5.Text += $"y2= [{minimumY3}; {maximumY3}]\n";
            richTextBox5.Text += $"y3= [{minimumY4}; {maximumY4}]\n\n";

            richTextBox5.Text += $"Limits for variable Bo\n";
            richTextBox5.Text += $"y0= [{minimumY1 / 2.4241}; {maximumY1 * 1.12}]\n";
            richTextBox5.Text += $"y1= [{minimumY2 * -2.53}; {maximumY2 * 23}]\n";
            richTextBox5.Text += $"y2= [{minimumY3 * -0.331}; {maximumY3 * 12}]\n";
            richTextBox5.Text += $"y3= [{minimumY4 * 0.4412}; {maximumY4 * 3.4}]\n\n";

            var charts = new Chart[] { chart1, chart2, chart3, chart4 };
            var errors = new TextBox[] { textBox9, textBox10, textBox11, textBox12 };
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
                errors[m].Text = max.ToString();
            }
        }       
        public void ReadFromFile(string file, double[,] x1, double[,] x2, double[,] x3, double[,] y)
        {
            string[] s = File.ReadAllLines(file);
            if (s.Length != x1.GetLength(0) + 1) 
                throw new Exception("Не співпадають розмірності !");

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
    
        public double T(double x, int m)
        {
            if (m == 0)
                return 1;
            if (m == 1)
                return x;
            return 2 * x * T(x, m - 1) - T(x, m - 2);
        }
        public double Tz(double x, int m)
        {
            if (m == 0)
                return 0.5;
            return T(2 * x - 1, m);
        }

        public void Normalization(double[,] matrix)
        {
            int size1 = matrix.GetLength(0);
            int size2 = matrix.GetLength(1);

            for (int i = 0; i < size2; ++i)
            {
                double min = matrix[0, i];
                double max = matrix[0, i];
                for (int j = 0; j < size1; ++j)
                {
                    if (matrix[j, i] > max)
                        max = matrix[j, i];
                    if (matrix[j, i] < min)
                        min = matrix[j, i];
                }
                for (int j = 0; j < size1; ++j)
                {
                    matrix[j, i] = (matrix[j, i] - min) / (max - min);
                }
            }
        }

        public double[] b(double[,] y)
        {
            int size1 = y.GetLength(0);
            int size2 = y.GetLength(1);

            double[] res = new double[size1];
            for (int i = 0; i < size1; ++i)
            {
                double min = y[i, 0];
                double max = y[i, 0];
                for (int j = 0; j < size2; ++j)
                {
                    if (y[i, j] > max)
                        max = y[i, j];
                    if (y[i, j] < min)
                        min = y[i, j];
                }
                res[i] = (max + min) / 2;
            }
            return res;
        }

        private void openToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (openFileDialog1.ShowDialog() == DialogResult.Cancel)
                return;
            string file = openFileDialog1.FileName;
            
            n1 = Int32.Parse(textBox1.Text);
            n2 = Int32.Parse(textBox2.Text);
            n3 = Int32.Parse(textBox3.Text);
            n = Int32.Parse(textBox4.Text);
            size = Int32.Parse(textBox8.Text);

            x1 = new double[size, n1];
            x2 = new double[size, n2];
            x3 = new double[size, n3];
            y = new double[size, n];

            ReadFromFile(file, x1, x2, x3, y);
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        public double[] Lambda(int p1, int p2, int p3)
        {
            double[] matrix_b = b(y);

            //int m = n1 * (p1 + 1) + n2 * (p2 + 1) + n3 * (p3 + 1);
            int m = n1 * p1 + n2 * p2 + n3 * p3 + 1;
            Matrix<double> matrix = Matrix<double>.Build.Dense(size, m);
            for (int i = 0; i < size; ++i) 
            {
                matrix[i, 0] = Math.Log(1.5) ;// 0.5 * (n1 + n2 + n3)
                for (int j = 0; j < n1; ++j)
                {
                    for (int p = 1; p < p1 + 1; ++p)
                    {
                        matrix[i, j * p1 + p] = Math.Log(1e-9+1 + Tz(x1[i, j], p));
                    }
                }
                for (int j = 0; j < n2; ++j)
                {
                    for (int p = 1; p < p2 + 1; ++p)
                    {
                        matrix[i, n1 * p1 + j * p2 + p] = Math.Log(1e-9 + 1 + Tz(x2[i, j], p));
                    }
                }
                for (int j = 0; j < n3; ++j)
                {
                    for (int p = 1; p < p3 + 1; ++p)
                    {
                        matrix[i, n1 * p1 + n2 * p2 + j * p3 + p] = Math.Log(1e-9 + 1 + Tz(x3[i, j], p));
                    }
                }
            }

            Matrix<double> vector = Matrix<double>.Build.Dense(size, 1);
            for (int i = 0; i < size; ++i)
            {
                vector[i, 0] = Math.Log(1e-9 + 1 + matrix_b[i]);
            }
            
            Matrix<double> resVector = (matrix.Transpose() * matrix).Inverse() * (matrix.Transpose() * vector);
           
            double[] res = new double[m];
            for (int i = 0; i < m; ++i)
            {
                res[i] = resVector[i, 0];
            }

            //Matrix<double> m1 = matrix.Transpose() * matrix;
            

            //using (StreamWriter sw = new StreamWriter("A1.txt"))
            //{
            //    //sw.WriteLine(m1.Determinant().ToString());
            //    for (int i = 0; i < m1.RowCount; ++i)
            //    {
            //        for (int j = 0; j < m1.ColumnCount; ++j)
            //        {
            //            sw.Write(String.Format("{0}\t", m1[i, j]));
            //        }
            //        sw.WriteLine();
            //    }
            //}

            //using (StreamWriter sw = new StreamWriter("b.txt"))
            //{
            //    for (int i = 0; i < size; ++i)
            //    {
            //        sw.Write(String.Format("{0}\t", matrix_b[i]));
            //        sw.WriteLine();
            //    }
            //}

            return res;
        }
        
        public double[,] Psi(double[,] Lambda, double[,] x)
        {
            int m = x.GetLength(1);
            double[,] res = new double[size, m];
            for (int i = 0; i < size; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    double sum = 0;
                    for (int p = 0; p < Lambda.GetLength(1); ++p)
                    {
                        sum += Lambda[j, p] *Math.Log(1e-9 + 1 +Tz(x[i, j], p));
                    }
                    res[i, j] = Math.Exp(sum)-1;
                }
            }
            return res;
        }
        public double[,] A(double[,] psi, double[,] y)
        {
            Matrix<double> _psi = Matrix<double>.Build.DenseOfArray(psi);
            Matrix<double> _y = Matrix<double>.Build.DenseOfArray(y);
            for(int i = 0; i<y.GetLength(0); ++i)
            {
                for (int j = 0; j<y.GetLength(1); ++j)
                {
                    _y[i,j] = Math.Log(1+ 1e-9 + _y[i,j]);
                }
            }
            for (int i = 0; i < psi.GetLength(0); ++i)
            {
                for (int j = 0; j < psi.GetLength(1); ++j)
                {
                    _psi[i, j] = Math.Log(1e-9 + 1 +_psi[i, j]);
                }
            }


            Matrix<double> resA = (_psi.Transpose() * _psi).Inverse() * (_psi.Transpose() * _y);
            double[,] res = new double[y.GetLength(1), psi.GetLength(1)];

            for (int i = 0; i < y.GetLength(1); ++i) 
            {
                for (int j = 0; j < psi.GetLength(1); ++j)
                {
                    res[i, j] = resA[j, i];
                }
            }
            return res;
        }

        public double[] Fi(double[,] A, double[,] Psi, int index)
        {
            double[] res = new double[size];
            for (int i = 0; i < size; ++i) 
            {
                double sum = 0;
                for (int j = 0; j < A.GetLength(1); ++j)
                {
                    sum += A[index, j] * Math.Log(Psi[i, j]+1);
                }
                res[i] = Math.Exp(sum)-1;
            }
            return res;
        }
        public double[] C(double[,] A1, double[,] A2, double[,] A3, double[,] Psi1, double[,] Psi2, double[,] Psi3, double[,] Y, int index)
        {
            Matrix<double> _Y = Matrix<double>.Build.Dense(size, 1);
            for (int i = 0; i < size; ++i)
            {
                _Y[i, 0] = Math.Log(Y[i, index]+1);
            }
            Matrix<double> matrix = Matrix<double>.Build.Dense(size, 3);
            for (int i = 0; i < size; ++i)
            {
                matrix[i, 0] = Math.Log(1 + 1e-9 + Fi(A1, Psi1, index)[i]);
                matrix[i, 1] = Math.Log(1 + 1e-9 + Fi(A2, Psi2, index)[i]);
                matrix[i, 2] = Math.Log(1 + 1e-9 + Fi(A3, Psi3, index)[i]);
            }

            Matrix<double> resC = (matrix.Transpose() * matrix).Inverse() * (matrix.Transpose() * _Y);

            double[] res = new double[3];
            for (int i = 0; i < 3; ++i)
            {
                res[i] = resC[i, 0];
            }
            return res;
        }
    }
}
