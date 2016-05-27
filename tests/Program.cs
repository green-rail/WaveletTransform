using System;
using System.Collections.Generic;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using static SignalGenerator.SignalGenerator;
using WaveletTransform;
using System.Drawing;

namespace Tests
{
    class Program
    {
        static void Main(string[] args)
        {
            Program prog = new Program();
            List<Wavelet> Daublets = WaveletConstructor.CreateAllDaubechies();
            List<Wavelet> Symlets = WaveletConstructor.CreateAllSymlets();
            //prog.Test1();
            prog.TestGetSet2d(Symlets);
            //prog.Test2D(Symlets);
            //prog.TestGetSet(Daublets);
            
            /*
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            TimeSpan ts;

            sw.Start();
            prog.TestSafe(Symlets);
            ts = sw.Elapsed;
            string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            Console.WriteLine("Test Safe finished in " + elapsedTime);
            Console.WriteLine();

            sw.Restart();
            prog.TestFast(Symlets);
            ts = sw.Elapsed;
            elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
            Console.WriteLine("Test Fast finished in " + elapsedTime);
            */
            Console.ReadKey();
        }

        void Test1()
        {
            string filePath = "test1.csv";
            string delimiter = ",";
            StringBuilder sb = new StringBuilder();
            Double[] TestSignalOne = CreateSignal(SignalList.Doppler, 0, 1024);
            Double[] TestSignalTwo = new Double[TestSignalOne.Length];
            List<Wavelet> Daublets = WaveletConstructor.CreateAllDaubechies();
            Transform.Forward1D(TestSignalOne, out TestSignalTwo, Daublets[1], 1);
            Transform.Inverse1D(TestSignalTwo, out TestSignalOne, Daublets[1], 1);
            for (int i= 0; i < TestSignalOne.Length; i++)
            {
                sb.AppendLine(string.Join(delimiter, TestSignalOne[i]));
            }
            File.WriteAllText(filePath, sb.ToString());
        }
        void TestSafe(List<Wavelet> wavelets)
        {
            Int32 Level = 5;
            Int32 len = 2048;
            //Int32 len = 32768;
            //Int32 len = 65536;
            //Int32 len = 131072;
            Double[] signal;
            Double[] signal2;
            Double SumError1 = 0;
            Double SumError2 = 0;
            Double SumError3 = 0;
            signal = CreateSignal(SignalList.Doppler, 0, len);
            foreach (Wavelet wavelet in wavelets)
            {
                Transform.Forward1D(signal, out signal2, wavelet, Level);
                Transform.Inverse1D(signal2, out signal2, wavelet, Level);
                SumError1 = AuxFunctions.ErrorNorm(signal, signal2);
                Console.WriteLine(wavelet.Name + " " + SumError1);
            }
            signal = CreateSignal(SignalList.Mixed, 0, len);
            foreach (Wavelet wavelet in wavelets)
            {
                Transform.Forward1D(signal, out signal2, wavelet, Level);
                Transform.Inverse1D(signal2, out signal2, wavelet, Level);
                SumError2 = AuxFunctions.ErrorNorm(signal, signal2);
                Console.WriteLine(wavelet.Name + " " + SumError2);
            }
            signal = CreateSignal(SignalList.Rectangular, 0, len);
            foreach (Wavelet wavelet in wavelets)
            {
                Transform.Forward1D(signal, out signal2, wavelet, Level);
                Transform.Inverse1D(signal2, out signal2, wavelet, Level);
                SumError3 = AuxFunctions.ErrorNorm(signal, signal2);
                Console.WriteLine(wavelet.Name + " " + SumError3);
            }
        }

        void TestFast(List<Wavelet> wavelets)
        {
            Int32 Level = 5;
            Int32 len = 2048;
            //Int32 len = 32768;
            //Int32 len = 65536;
            //Int32 len = 131072;
            Double[] signal;
            Double[] signal2;
            Double SumError1 = 0;
            Double SumError2 = 0;
            Double SumError3 = 0;
            signal = CreateSignal(SignalList.Doppler, 0, len);
            foreach (Wavelet wavelet in wavelets)
            {
                Transform.FastForward1d(signal, out signal2, wavelet, Level);
                Transform.FastInverse1d(signal2, out signal2, wavelet, Level);
                SumError1 = AuxFunctions.ErrorNorm(signal, signal2);
                Console.WriteLine(wavelet.Name + " " + SumError1);
            }
            signal = CreateSignal(SignalList.Mixed, 0, len);
            foreach (Wavelet wavelet in wavelets)
            {
                Transform.FastForward1d(signal, out signal2, wavelet, Level);
                Transform.FastInverse1d(signal2, out signal2, wavelet, Level);
                SumError2 = AuxFunctions.ErrorNorm(signal, signal2);
                Console.WriteLine(wavelet.Name + " " + SumError2);
            }
            signal = CreateSignal(SignalList.Rectangular, 0, len);
            foreach (Wavelet wavelet in wavelets)
            {
                Transform.FastForward1d(signal, out signal2, wavelet, Level);
                Transform.FastInverse1d(signal2, out signal2, wavelet, Level);
                SumError3 = AuxFunctions.ErrorNorm(signal, signal2);
                Console.WriteLine(wavelet.Name + " " + SumError3);
            }
        }

        void Test2D(List<Wavelet> wavelets)
        {
            Bitmap OrImage = new Bitmap("LENA1.bmp");
            Double[,] ImageVals = new Double[OrImage.Width, OrImage.Height];
            Color Col;

            for (int i = 0; i < OrImage.Width; i++)
            {
                for (int j = 0; j < OrImage.Height; j++)
                {
                    Col = OrImage.GetPixel(i, j);
                    ImageVals[i, j] = AuxFunctions.Scale(0, 255, -1, 1, Col.R);
                }
            }
            Transform.FastForward2d(ImageVals, out ImageVals, wavelets[1], 3);
            Transform.FastInverse2d(ImageVals, out ImageVals, wavelets[1], 3);
            for (int i = 0; i < OrImage.Width; i++)
            {
                for (int j = 0; j < OrImage.Height; j++)
                {
                    int buff = (Int32)AuxFunctions.Scale(-1, 1, 0, 255, ImageVals[i,j]);
                    Col = Color.FromArgb(buff, buff, buff);
                    OrImage.SetPixel(i, j, Col);
                }
            }

            OrImage.Save("processed.bmp");
        }

        public void TestGetSet(List<Wavelet> waves)
        {
            Double[] Ssig = { 6, 3, 0, 8 };
            Transform.FastForward1d(Ssig, out Ssig, waves[0], 2);
            Double[] Details1 = Transform.GetDetailOfLevel(Ssig, 2, 1);
            Double[] Details2 = Transform.GetDetailOfLevel(Ssig, 2, 2);
            foreach (Double val in Ssig) Console.WriteLine(val);
            Console.WriteLine();
            foreach (Double val in Details1) Console.WriteLine(val);
            foreach (Double val in Details2) Console.WriteLine(val);
            Transform.SetDetailOfLevel(Ssig, new Double[] { 0.666, 0.777});
            Console.WriteLine("We just set up new details of level one...");
            Details1 = Transform.GetDetailOfLevel(Ssig, 2, 1);
            foreach (Double val in Details1) Console.WriteLine(val);
            Console.WriteLine("and here's new coefs:");
            foreach (Double val in Ssig) Console.WriteLine(val);
        }
        public void TestGetSet2d(List<Wavelet> wavelets)
        {
            Double[,] Signal = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 }, { 8, 9, 10, 11 }, { 12, 13, 14, 15 } };
            Double[] details = Transform.GetAllDetail(Signal, 1);
            for (int i = 0; i < details.Length; i++)
            {
                details[i] = 666;
            }
            Transform.SetAllDetail(details, Signal, 1);

            /*
            foreach(Double det in details)
            {
                Console.Write(det + " ");
            }

            Console.ReadKey();
            */
            foreach (Double det in Signal)
            {
                Console.WriteLine(det + " ");
            }

        }
    }
}
