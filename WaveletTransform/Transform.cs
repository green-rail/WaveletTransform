using System;

namespace WaveletTransform
{
    public static class Transform
    {
        public static unsafe void FastForward1d(Double[] input, out Double[] output, Wavelet wavelet, Int32 Level)
        {
            Int32 Len = wavelet.DecompositionLow.Length;
            Int32 CircleInd;
            output = new Double[input.Length];
            Double[] Buff = new Double[input.Length];
            Buffer.BlockCopy(input, 0, output, 0, input.Length * 8);
            Double BufScal = 0;
            Double BufDet = 0;
            Double* DecLow = stackalloc Double[Len];
            Double* DecHigh = stackalloc Double[Len];
            for (int i = 0; i < Len; i++)
            {
                DecLow[i] = wavelet.DecompositionLow[i];
                DecHigh[i] = wavelet.DecompositionHigh[i];
            }

            fixed (Double* pout = output, pbuf = Buff){

                for (int level = 0; level < Level; level++)
                {
                    Int32 Bound = input.Length >> level;
                    Int32 StartIndex = -((Len >> 1) - 1);
                    Buffer.BlockCopy(output, 0, Buff, 0, Bound * 8);

                    for (int i = 0; i < Bound >> 1; i++)
                    {
                        for (int j = StartIndex, k = 0; k < Len; j++, k++)
                        {
                            if ((StartIndex < 0) || j >= Bound) CircleInd = ((j % Bound) + Bound) % Bound;
                            else CircleInd = j;
                            BufScal += DecLow[k] * pout[CircleInd];
                            BufDet += DecHigh[k] * pout[CircleInd];
                        }
                        StartIndex += 2;
                        pbuf[i] = BufScal;
                        pbuf[i + (Bound >> 1)] = BufDet;
                        BufScal = 0;
                        BufDet = 0;
                    }
                    Buffer.BlockCopy(Buff, 0, output, 0, Bound * 8);
                }
            }
        }

        public static unsafe void FastInverse1d(Double[] input, out Double[] output, Wavelet wavelet, Int32 startLevel)
        {
            Int32 Len = wavelet.ReconstructionLow.Length;
            Int32 CircleInd;
            output = new Double[input.Length];
            Double[] Buff = new Double[input.Length];
            Double[] BufferLow = new Double[input.Length];
            Double[] BufferHigh = new Double[input.Length];
            Buffer.BlockCopy(input, 0, output, 0, input.Length * 8);
            Double Buf = 0;
            Double* RecLow = stackalloc Double[Len];
            Double* RecHigh = stackalloc Double[Len];
            for (int i = 0; i < Len; i++)
            {
                RecLow[i] = wavelet.ReconstructionLow[i];
                RecHigh[i] = wavelet.ReconstructionHigh[i];
            }

            fixed (Double* pbuf = Buff, pLow = BufferLow, pHigh = BufferHigh){

                for (int level = startLevel; level > 0; level--)
                {
                    Int32 Bound = input.Length >> level;
                    Int32 StartIndex = -((Len >> 1) - 1);
                    for (int i = 0, j = 0; i < Bound << 1; i += 2, j++)
                    {
                        pLow[i] = 0;
                        pHigh[i] = 0;
                        pLow[i + 1] = output[j];
                        pHigh[i + 1] = output[Bound + j];
                    }

                    for (int i = 0; i < Bound << 1; i++)
                    {
                        for (int j = StartIndex, k = 0; k < Len; j++, k++)
                        {
                            if ((StartIndex < 0) || j >= (Bound << 1)) CircleInd = (j % (Bound << 1) + (Bound << 1)) % (Bound << 1);
                            else CircleInd = j;
                            Buf += RecLow[k] * pLow[CircleInd] + RecHigh[k] * pHigh[CircleInd];
                        }
                        StartIndex += 1;
                        pbuf[i] = Buf;
                        Buf = 0;
                    }
                    Buffer.BlockCopy(Buff, 0, output, 0, Bound * 16);
                }
            }
        }


        private static Double[] StepForward(Double[] input, Wavelet wavelet, int CurrentLevel)
        {
            Int32 Bound = input.Length >> CurrentLevel;
            Int32 Len = wavelet.DecompositionLow.Length;
            Int32 StartIndex = -((Len >> 1) - 1);
            Double[] Output = new Double[input.Length];
            Array.Copy(input, Bound, Output, Bound, input.Length - Bound);
            Double BufScal = 0;
            Double BufDet = 0;
            Int32 CircleInd;
            for (int i = 0, r = 0; i < Bound; i+=2, r++)
            {
                for (int j = StartIndex, k = 0; k < Len; j++, k++)
                {
                    if ((StartIndex < 0) || j >= Bound) CircleInd = ((j % Bound) + Bound) % Bound;
                    else CircleInd = j;
                    BufScal += wavelet.DecompositionLow[k] * input[CircleInd];
                    BufDet += wavelet.DecompositionHigh[k] * input[CircleInd];
                }
                StartIndex += 2;
                Output[r] = BufScal;
                Output[r + (Bound >> 1)] = BufDet;
                BufScal = 0;
                BufDet = 0;
            }
            return Output;
        }
        public static void Forward1D(Double[] input, out Double[] output, Wavelet wavelet, Int32 Level)
        {
            bool f = (input.Length & (input.Length - 1)) == 0; 
            if (!f || input.Length == 0) throw new Exception("input length must be a power of two");

            output = new Double[input.Length];
            Array.Copy(input, output, input.Length);
            for (int i = 0; i < Level; i++)
            {
                output = StepForward(output, wavelet, i);
            }
        }
        private static Double[] StepInverse(Double[] input, Wavelet wavelet, int CurrentLevel)
        {
            Int32 Bound = input.Length >> CurrentLevel;
            Int32 Len = wavelet.ReconstructionLow.Length;
            Int32 StartIndex = -((Len >> 1) - 1);
            Double[] Output = new Double[input.Length];
            Double[] BuffLow = new Double[Bound << 1];
            Double[] BuffHi = new Double[Bound << 1];
            Array.Copy(input, 0, Output, 0, input.Length);
            Double BufScal = 0;
            Int32 CircleInd;
            for (int i = 0, j = 0; i < Bound << 1; i+=2, j++)
            {
                BuffLow[i] = 0;
                BuffHi[i] = 0;
                BuffLow[i + 1] = input[j];
                BuffHi[i + 1] = input[Bound + j];
            }
            for (int i = 0; i < Bound << 1; i++)
            {
                for (int j = StartIndex, k = 0; k < Len; j++, k++)
                {
                    if ((StartIndex < 0) || j >= (Bound << 1)) CircleInd = (j % (Bound << 1) + (Bound << 1)) % (Bound << 1);
                    else CircleInd = j;
                    BufScal += wavelet.ReconstructionLow[k] * BuffLow[CircleInd] + wavelet.ReconstructionHigh[k] * BuffHi[CircleInd];
                }
                StartIndex += 1;
                Output[i] = BufScal;
                BufScal = 0;
            }
            return Output;
        }
        public static void Inverse1D(Double[] input, out Double[] output, Wavelet wavelet, Int32 StartLevel)
        {
            bool f = (input.Length & (input.Length - 1)) == 0; 
            if (!f || input.Length == 0) throw new Exception("input length must be a power of two");

            output = new Double[input.Length];
            Array.Copy(input, output, input.Length);
            for (int i = StartLevel; i > 0; i--)
            {
                output = StepInverse(output, wavelet, i);
            }
        }


        private static unsafe void FastStepForward(ref Double[] input, Double[] buff, Int32 len, Int32 CurrentLevel, Double* DecLow, Double* DecHigh)
        {
            Int32 CircleInd;
            Double BufScal = 0;
            Double BufDet = 0;
            Int32 Bound = input.Length >> CurrentLevel;
            Int32 StartIndex = -((len >> 1) - 1);
            fixed (Double* pinp = input, pbuf = buff)
            {
                for (int i = 0; i < Bound >> 1; i++)
                {
                    for (int j = StartIndex, k = 0; k < len; j++, k++)
                    {
                        if ((StartIndex < 0) || j >= Bound) CircleInd = ((j % Bound) + Bound) % Bound;
                        else CircleInd = j;
                        BufScal += DecLow[k] * pinp[CircleInd];
                        BufDet += DecHigh[k] * pinp[CircleInd];
                    }
                    StartIndex += 2;
                    pbuf[i] = BufScal;
                    pbuf[i + (Bound >> 1)] = BufDet;
                    BufScal = 0;
                    BufDet = 0;
                }
                Buffer.BlockCopy(buff, 0, input, 0, Bound * 8);
            }
        }

        private static unsafe void FastStepInverse(ref Double[] input, Double[] buffLow, Double[] buffHigh, Int32 len, Int32 CurrentLevel, Double* RecLow, Double* RecHigh)
        {
            Int32 Bound = input.Length >> CurrentLevel;
            Int32 StartIndex = -((len >> 1) - 1);
            Double Buf = 0;
            Int32 CircleInd;
            fixed (Double* pinp = input, pbufLow = buffLow, pbufHi = buffHigh)
            {
                for (int i = 0, j = 0; i < Bound << 1; i+=2, j++)
                {
                    pbufLow[i] = 0;
                    pbufHi[i] = 0;
                    pbufLow[i + 1] = pinp[j];
                    pbufHi[i + 1] = pinp[Bound + j];
                }
                for (int i = 0; i < Bound << 1; i++)
                {
                    for (int j = StartIndex, k = 0; k < len; j++, k++)
                    {
                        if ((StartIndex < 0) || j >= (Bound << 1)) CircleInd = (j % (Bound << 1) + (Bound << 1)) % (Bound << 1);
                        else CircleInd = j;
                        Buf += RecLow[k] * pbufLow[CircleInd] + RecHigh[k] * pbufHi[CircleInd];
                    }
                    StartIndex += 1;
                    pinp[i] = Buf;
                    Buf = 0;
                }
            }
        }

        public static unsafe void FastForward2d (Double[,] input, out Double[,] output, Wavelet wavelet, Int32 Level)
        {
            Int32 DataLen = input.GetLength(0);
            Int32 Len = wavelet.DecompositionHigh.Length;
            Int32 Bound;
            output = new Double[DataLen, DataLen];
            Double[] buff = new Double[DataLen];
            Double[] buffData = new Double[DataLen];

            Double* DecLow = stackalloc Double[Len];
            Double* DecHigh = stackalloc Double[Len];
            for (int i = 0; i < Len; i++)
            {
                DecLow[i] = wavelet.DecompositionLow[i];
                DecHigh[i] = wavelet.DecompositionHigh[i];
            }

            for (int i = 0; i < DataLen; i++)
            {
                for (int j = 0; j < DataLen; j++)
                {
                    output[i, j] = input[i, j];
                }
            }

            for (int lev = 0; lev < Level; lev++)
            {
                Bound = DataLen >> lev;
                for (int i = 0; i < Bound; i++)
                {
                    for (int j = 0; j < Bound; j++) buffData[j] = output[i, j];
                    FastStepForward(ref buffData, buff, Len, lev, DecLow, DecHigh);
                    for (int j = 0; j < Bound; j++) output[i, j] = buffData[j]; 
                }
                for (int j = 0; j < Bound; j++)
                {
                    for (int i = 0; i < Bound; i++) buffData[i] = output[i, j];
                    FastStepForward(ref buffData, buff, Len, lev, DecLow, DecHigh);
                    for (int i = 0; i < Bound; i++) output[i, j] = buffData[i];
                }
            }
        }

        public static unsafe void FastInverse2d (Double[,] input, out Double[,] output, Wavelet wavelet, Int32 Level)
        {
            Int32 DataLen = input.GetLength(0);
            Int32 Len = wavelet.DecompositionHigh.Length;
            Int32 Bound;
            output = new Double[DataLen, DataLen];
            Double[] buffData = new Double[DataLen];
            Double[] buffLow = new Double[DataLen];
            Double[] buffHigh = new Double[DataLen];

            Double* RecLow = stackalloc Double[Len];
            Double* RecHigh = stackalloc Double[Len];
            for (int i = 0; i < Len; i++)
            {
                RecLow[i] = wavelet.ReconstructionLow[i];
                RecHigh[i] = wavelet.ReconstructionHigh[i];
            }

            for (int i = 0; i < DataLen; i++)
            {
                for (int j = 0; j < DataLen; j++)
                {
                    output[i, j] = input[i, j];
                }
            }

            for (int lev = Level; lev > 0; lev--)
            {
                Bound = DataLen >> lev;
                for (int j = 0; j < Bound << 1; j++)
                {
                    for (int i = 0; i < Bound << 1; i++) buffData[i] = output[i, j];
                    FastStepInverse(ref buffData, buffLow, buffHigh, Len, lev, RecLow, RecHigh);
                    for (int i = 0; i < Bound << 1; i++) output[i, j] = buffData[i]; 
                }
                for (int i = 0; i < Bound << 1; i++)
                {
                    for (int j = 0; j < Bound << 1; j++) buffData[j] = output[i, j];
                    FastStepInverse(ref buffData, buffLow, buffHigh, Len, lev, RecLow, RecHigh);
                    for (int j = 0; j < Bound << 1; j++) output[i, j] = buffData[j];
                }
            }
        }

        public static Double[] GetAllDetail(Double[] input, Int32 currentLevel)
        {
            Int32 Shift = input.Length >> currentLevel;
            Double[] output = new Double[input.Length - Shift];
            Array.Copy(input, Shift, output, 0, input.Length - Shift);
            return output;
        }
        public static Double[] GetAllDetail(Double[,] input, Int32 currentLevel)
        {
            Int32 Shift = input.Length >> currentLevel;
            Int32 Len = input.GetLength(0) >> currentLevel;
            Double[] output = new Double[input.GetLength(0)* input.GetLength(0) - (Len * Len) ];
            Int32 OutIndex = 0;

            for (int lev = currentLevel; lev > 0; lev--)
            {
                Int32 Bound = input.GetLength(0) >> lev - 1;
                Len = input.GetLength(0) >> lev;
                for (int i = 0; i < Len; i++)
                {
                    for (int j = Len; j < Bound; j++)
                    {
                        output[OutIndex] = input[i, j];
                        OutIndex++;
                    }
                }
                for (int i = Len; i < Bound; i++)
                {
                    for (int j = 0; j < Len; j++)
                    {
                        output[OutIndex] = input[i, j];
                        OutIndex++;
                    }
                }
                for (int i = Len; i < Bound; i++)
                {
                    for (int j = Len; j < Bound; j++)
                    {
                        output[OutIndex] = input[i, j];
                        OutIndex++;
                    }
                }
            }
            return output;
        }

        public static void SetAllDetail(Double[] input, Double[,] coefs, Int32 currentLevel)
        {
            Int32 Shift = coefs.Length >> currentLevel;
            Int32 Len = coefs.GetLength(0) >> currentLevel;
            Int32 InpIndex = 0;

            //for (int lev = 1; lev > currentLevel; lev++)
            for (int lev = currentLevel; lev > 0; lev--)
            {
                Int32 Bound = coefs.GetLength(0) >> lev - 1;
                Len = coefs.GetLength(0) >> lev;
                for (int i = 0; i < Len; i++)
                {
                    for (int j = Len; j < Bound; j++)
                    {
                        coefs[i, j] = input[InpIndex]; 
                        InpIndex++;
                    }
                }
                for (int i = Len; i < Bound; i++)
                {
                    for (int j = 0; j < Len; j++)
                    {
                        coefs[i, j] = input[InpIndex]; 
                        InpIndex++;
                    }
                }
                for (int i = Len; i < Bound; i++)
                {
                    for (int j = Len; j < Bound; j++)
                    {
                        coefs[i, j] = input[InpIndex]; 
                        InpIndex++;
                    }
                }
            }
        }

        public static Double[] GetDetailOfLevel(Double[] input, Int32 currentLevel, Int32 level)
        {
            Int32 Len = input.Length >> level;
            Double[] output = new Double[Len];
            Array.Copy(input, Len, output, 0, Len);
            return output;
        }
        public static Double[] GetScaling(Double[] input, Int32 currentLevel)
        {
            Int32 Len = input.Length >> currentLevel;
            Double[] output = new Double[Len];
            Array.Copy(input, output, Len);
            return output;
        }
        public static void SetDetailOfLevel(Double[] input, Double[] inpDetails)
        {
            Array.Copy(inpDetails, 0, input, inpDetails.Length, inpDetails.Length);
        }
    }
}
