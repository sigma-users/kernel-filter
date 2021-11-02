using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Runtime.InteropServices;
using System.Collections;
using System.Threading;
using System.ComponentModel;
using System.Data;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;

namespace KernelFilter
{
	unsafe class Program
	{
		public static bool ReadBinaryFile(string path, int hpc, int vpc, out double[,] data)
		{
			data = null;
			try
			{
				data = new double[vpc, hpc];
				BinaryReader br = new BinaryReader(new FileStream(path, FileMode.Open));
				for (int i = 0; i < vpc; i++)
				{ for (int j = 0; j < hpc; j++) { data[i, j] = br.ReadDouble(); } }
				br.Close();
				return true;
			}
			catch { return false; }
		}

		public static void CreateScalarBmp(double[,] image, double cnt, out Bitmap bmp)
		{
			int height = image.GetLength(0);
			int width = image.GetLength(1);

			double min, max;
			min = image[0, 0];
			max = min;
			for(int i = 0; i < height; i++)
			{
				for(int j = 0; j < width; j++)
				{
					if(min > image[i, j]) { min = image[i, j]; }
					if(max < image[i, j]) { max = image[i, j]; }
				}
			}

			byte[] data = new byte[4 * height * width];
			fixed (byte* pdata = data)
			{
				fixed (double* pimg = image)
				{
					FImgScalarPlotData(0, width, height, pimg, 
						min, (max - min) / cnt + min, pdata);
				}
			}

			bmp = new Bitmap(width, height, PixelFormat.Format32bppRgb);
			Rectangle rect = new Rectangle(0, 0, bmp.Width, bmp.Height);
			BitmapData bdata = bmp.LockBits(rect, ImageLockMode.ReadWrite, bmp.PixelFormat);
			IntPtr ptr = bdata.Scan0;
			Marshal.Copy(data, 0, ptr, data.Length);
			bmp.UnlockBits(bdata);
		}

		public static void CreateVectorBmp(double[,] vx, double[,] vy, double cnt, out Bitmap bmp)
		{
			int height = vx.GetLength(0);
			int width = vx.GetLength(1);

			double min, max;
			min = vx[0, 0] * vx[0, 0] + vy[0, 0] * vy[0, 0];
			max = min;
			for (int i = 0; i < height; i++)
			{
				for (int j = 0; j < width; j++)
				{
					double v = vx[i, j] * vx[i, j] + vy[i, j] * vy[i, j];
					if (min > v) { min = v; }
					if (max < v) { max = v; }
				}
			}
			min = Math.Sqrt(min);
			max = Math.Sqrt(max);

			byte[] data = new byte[4 * height * width];
			fixed (byte* pdata = data)
			{
				fixed (double* pvx = vx, pvy = vy)
				{
					FImgVectorPlotData(width, height, pvx, pvy, 
						0.0, (max - min) / cnt, pdata);
				}
			}

			bmp = new Bitmap(width, height, PixelFormat.Format32bppRgb);
			Rectangle rect = new Rectangle(0, 0, bmp.Width, bmp.Height);
			BitmapData bdata = bmp.LockBits(rect, ImageLockMode.ReadWrite, bmp.PixelFormat);
			IntPtr ptr = bdata.Scan0;
			Marshal.Copy(data, 0, ptr, data.Length);
			bmp.UnlockBits(bdata);
		}

		static unsafe void CalcFilteredImage(
			double[] p, double[] q, double[,] org, out double[,] dest)
		{
			int pave = 1;
			double[] have = new double[pave];
			double[] vave = new double[pave];
			double[] wave = new double[pave];
			wave[0] = 1.0;

			int hpc = org.GetLength(1);
			int vpc = org.GetLength(0);
			dest = new double[vpc, hpc];
			int[,] mask = new int[vpc, hpc];
			fixed (int* pmask = mask)
			{
				fixed (double* pp = p, pq = q, porg = org, pdest = dest,
					phave = have, pvave = vave, pwave = wave)
				{
					FDApplyKernelFilterMini(
						0, hpc, vpc, pp, pq, porg, 5.0,
						pave, phave, pvave, pwave,
						pdest, pmask);
				}
			}
		}

		static void Main(string[] args)
		{
			//=======================================================================
			//  Main data processing
			//=======================================================================
			double[,] mainAdf, mainDx, mainDy;
			ReadBinaryFile("Main_DF_942_942.bin", 942, 942, out mainAdf);
			ReadBinaryFile("Main_Dx_942_942.bin", 942, 942, out mainDx);
			ReadBinaryFile("Main_Dy_942_942.bin", 942, 942, out mainDy);

			double[] mainP = new double[2];
			double[] mainQ = new double[2];
			mainP[0] = 23.0120177886009;
			mainP[1] = -30.4330035074283;
			mainQ[0] = -23.9449357015599;
			mainQ[1] = -14.6715391044128;

			double[,] mainFilteredDx, mainFilteredDy;
			CalcFilteredImage(mainP, mainQ, mainDx, out mainFilteredDx);
			CalcFilteredImage(mainP, mainQ, mainDy, out mainFilteredDy);

			Bitmap mainAdfBmp, mainMapBmp, mainFilteredMapBmp; ;
			CreateScalarBmp(mainAdf, 1.0, out mainAdfBmp);
			CreateVectorBmp(mainDx, mainDy, 1.0, out mainMapBmp);
			CreateVectorBmp(mainFilteredDx, mainFilteredDy, 2.0, out mainFilteredMapBmp);

			mainAdfBmp.Save("MainAdf.bmp");
			mainMapBmp.Save("MainMap.bmp");
			mainFilteredMapBmp.Save("MainFiltered.bmp");

			//=======================================================================
			//  RT data processing
			//=======================================================================
			double[,] rtAdf, rtDx, rtDy;
			ReadBinaryFile("RT_DF_912_912.bin", 912, 912, out rtAdf);
			ReadBinaryFile("RT_Dx_912_912.bin", 912, 912, out rtDx);
			ReadBinaryFile("RT_Dy_912_912.bin", 912, 912, out rtDy);

			double[] rtP = new double[2];
			double[] rtQ = new double[2];
			rtP[0] = 6.56942044171111;
			rtP[1] = 37.5130328013731;
			rtQ[0] = 26.5855697245923;
			rtQ[1] = -7.61413595413948;

			double[,] rtFilteredDx, rtFilteredDy;
			CalcFilteredImage(rtP, rtQ, rtDx, out rtFilteredDx);
			CalcFilteredImage(rtP, rtQ, rtDy, out rtFilteredDy);

			Bitmap rtAdfBmp, rtMapBmp, rtFilteredMapBmp; ;
			CreateScalarBmp(rtAdf, 1.0, out rtAdfBmp);
			CreateVectorBmp(rtDx, rtDy, 1.0, out rtMapBmp);
			CreateVectorBmp(rtFilteredDx, rtFilteredDy, 2.0, out rtFilteredMapBmp);

			rtAdfBmp.Save("RtAdf.bmp");
			rtMapBmp.Save("RtMap.bmp");
			rtFilteredMapBmp.Save("RtFiltered.bmp");

			//=======================================================================
			//  N2 data processing
			//=======================================================================
			double[,] n2Adf, n2Dx, n2Dy;
			ReadBinaryFile("N2_DF_860_860.bin", 860, 860, out n2Adf);
			ReadBinaryFile("N2_Dx_860_860.bin", 860, 860, out n2Dx);
			ReadBinaryFile("N2_Dy_860_860.bin", 860, 860, out n2Dy);

			double[] n2P = new double[2];
			double[] n2Q = new double[2];
			n2P[0] = 6.65511794208699;
			n2P[1] = 37.7653045871795;
			n2Q[0] = 26.5474447551317;
			n2Q[1] = -7.66953263523549;

			double[,] n2FilteredDx, n2FilteredDy;
			CalcFilteredImage(n2P, n2Q, n2Dx, out n2FilteredDx);
			CalcFilteredImage(n2P, n2Q, n2Dy, out n2FilteredDy);

			Bitmap n2AdfBmp, n2MapBmp, n2FilteredMapBmp; ;
			CreateScalarBmp(n2Adf, 1.0, out n2AdfBmp);
			CreateVectorBmp(n2Dx, n2Dy, 1.0, out n2MapBmp);
			CreateVectorBmp(n2FilteredDx, n2FilteredDy, 2.0, out n2FilteredMapBmp);

			n2AdfBmp.Save("N2Adf.bmp");
			n2MapBmp.Save("N2Map.bmp");
			n2FilteredMapBmp.Save("N2Filtered.bmp");
		}

		[DllImport("FKLib.dll")]
		public static extern unsafe void FDApplyKernelFilterMini(
			int filter_type, int hpc, int vpc, double* rp, double* rq,
			double* org, double edge, int pave, 
			double* have, double* vave, double* wave,
			double* image, int* mask);

		[DllImport("FKLib.dll")]
		public static extern unsafe void FImgScalarPlotData(
			int color_index, int width, int height, double* image,
			double smin, double smax, byte* gray);

		[DllImport("FKLib.dll")]
		public static extern unsafe void FImgVectorPlotData(
			int width, int height, double* vx, double* vy,
			double vmin, double vmax, byte* rgb);
	}
}
