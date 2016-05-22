using System;

namespace Tests
{
    public class AuxFunctions
    {
        public static double Median(double[] input) //Just gives median
        {
            double[] Buffer = new double[input.Length];
            for(int i = 0; i < input.Length; i++)
            {
                Buffer[i] = input[i];
            }

            Array.Sort(Buffer);
            if (Buffer.Length % 2 == 0)
            {
                return (Buffer[Buffer.Length / 2] + Buffer[Buffer.Length / 2 - 1]) / 2;
            }
            else return Buffer[Buffer.Length / 2];
        }

        public static double MedianOfModulus(double[] input)//Gives median of absolute values of the original vector
        {
            double[] Buffer = new double[input.Length];
            for (int i = 0; i < input.Length; i++)
            {
                Buffer[i] = Math.Abs(input[i]);
            }
            Array.Sort(Buffer);
            if (Buffer.Length % 2 == 0)
            {
                return (Buffer[Buffer.Length / 2] + Buffer[Buffer.Length / 2 - 1]) / 2;
            }
            else return Buffer[Buffer.Length / 2];
        }
        public static double MedianOfModulus(double[,] input)//Gives median of absolute values of the original vector доделать!
        {
            double[] Buffer = new double[input.Length];
            for (int i = 0; i < input.GetLength(0); i++)
            {
                for (int j = 0; j < input.GetLength(1); j++)
                {
                    Buffer[i] = Math.Abs(input[i, j]);
                }
            }
            Array.Sort(Buffer);
            if (Buffer.Length % 2 == 0)
            {
                return (Buffer[Buffer.Length / 2] + Buffer[Buffer.Length / 2 - 1]) / 2;
            }
            else return Buffer[Buffer.Length / 2];
        }

        public static double ErrorNorm(double[] clear, double[] filtrated) //Kinda like SNR, but opposite: ErrorNorm = Norm(Filtrated - Clear) / Norm(Clear) 
        {
            double[] ErrorArr = new double[clear.Length];
            for (int i = 0; i < clear.Length; i++)
            {
                ErrorArr[i] = filtrated[i] - clear[i];
            }
            return GetNorm(ErrorArr) / GetNorm(clear);
        }
        public static double ErrorNorm(double[,] clear, double[,] filtrated) //Kinda like SNR, but opposite: ErrorNorm = Norm(Filtrated - Clear) / Norm(Clear) 
        {
            double[,] ErrorArr = new double[clear.GetLength(0), clear.GetLength(1)];
            for (int i = 0; i < clear.GetLength(0); i++)
            {
                for(int j = 0; j < clear.GetLength(1); j++)
                {
                    ErrorArr[i,j] = filtrated[i,j] - clear[i,j];
                }
            }
            return GetNorm(ErrorArr) / GetNorm(clear);
        }

        private static double GetNorm(double[] input) //Euclidean norm
        {
            double Buff = 0;
            for (int i = 0; i < input.Length; i++)
            {
                Buff += input[i] * input[i];
            }
            return Math.Sqrt(Buff);
        }        
        private static double GetNorm(double[,] input) //Euclidean norm 2D
        {
            double Buff = 0;
            for (int i = 0; i < input.GetLength(0); i++)
            {
                for (int j = 0; j < input.GetLength(1); j++) Buff += input[i, j] * input[i, j];
            }
            return Math.Sqrt(Buff); 
        }        

        public static double Scale(double fromMin, double fromMax, double toMin, double toMax, double x)
        {
            if (fromMax - fromMin == 0) return 0;
            double value = (toMax - toMin) * (x - fromMin) / (fromMax - fromMin) + toMin;
            if (value > toMax)
            {
                value = toMax;
            }
            if (value < toMin)
            {
                value = toMin;
            }
            return value;
        }
    }
}
