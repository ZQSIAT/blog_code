//opencv版本:OpenCV2.7.13
//VS版本:VS2015
//Author:zqs
//time:2018-03-31 18:36:00
//name:LoG

#include<iostream>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>

#include<math.h>
#include<fstream>
#include<string>

using namespace std;
using namespace cv;

//两种卷积方式
//#define filter2D_convolution
#define my_convolution

//三种二值化方式
//#define Only_threshold
//#define threshold_neighborhood
#define neighborhood_negative_2
#define neighborhood_negative_1

//是否存储二值化结果
//#define store_result

//定义卷积核大小
#define N_3
//#define N_5
//#define N_9
//#define N_11

#ifdef N_3
#define Sigma 0.5
#define N 3
string dst_str_1 = "σ=0.5,N=3";
#endif

#ifdef N_5
#define Sigma 1.2
#define N 5
string dst_str_1 = "σ=1.2,N=5";
#endif

#ifdef N_9
#define Sigma 2.8
#define N 9
string dst_str_1 = "σ=2.8,N=9";
#endif

#ifdef N_11
#define Sigma 1.8
#define N 11
string dst_str_1 = "σ=1.8,N=11";
#endif

//存储图像为txt
void store_data(Mat &image, string txt_name)
{

	ofstream outfile;
	outfile.open(txt_name + ".txt"); //存放数据的文件名

	Mat input_image = image.clone();//复制实参到临时变量
	int n_r = input_image.rows; // 行数
	int n_c = input_image.cols * input_image.channels();//列数*通道数等于每一行元素的个数
	for (int j = 0; j<n_r; j++)
	{
		#ifdef filter2D_convolution
		uchar* data = input_image.ptr<uchar>(j);
		#else
		double* data = input_image.ptr<double>(j);
		#endif
		for (int i = 0; i<n_c; i++)
		{
			if (outfile.is_open())
			{

				outfile << (double)data[i] << " "; //程序中处理的数据
			}
			else
			{
				cout << "不能打开文件!" << endl;
			}
		}
		if (outfile.is_open())
		{
			outfile << endl;
		}
		else
		{
			cout << "不能打开文件!" << endl;
		}
	}
	outfile.close();

}
//得到LoG卷积模板，根据公式得出double型模板，然后调整整体位置，使模板内系数和为0。
//N为模板大小，其值为奇数。double sigma为高斯函数的标准差。
int getLoGModel(string str, Mat src, double sigma, int _threshold)
{
	if (N % 2 == 0)
	{
		cout << "Error: N is not a even!" << endl;
	}

	double model[N*N];
	double model_convolution[N][N];
	double sumLaplacian;

	double x, y, center_x, center_y;
	center_x = (N - 1) / 2;
	center_y = (N - 1) / 2;
	for (int j = 0; j<N; j++)
	{
		for (int i = 0; i<N; i++)
		{
			x = i - center_x;
			y = j - center_y;
			model[j*N + i] = -((x*x + y*y - 2 * sigma*sigma) / pow(sigma, 4))*(exp(-(x*x + y*y) / (2 * sigma*sigma)));
		}
	}
	double sum = 0.0;
	for (int i = 0; i<N*N; i++)
	{
		sum += model[i];
	}
	double average = sum / (double)(N*N);
	for (int i = 0; i<N*N; i++)
	{
		model[i] -= average;
	}
	//显示卷积模板1并赋值给model_convolution[N][N]
	cout << "N = " << endl;
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
		{
			model_convolution[i][j] = model[i*N + j];
			cout << model[i*N + j] << ' ';
		}
		cout << endl;
	}

	#ifdef filter2D_convolution
	Mat dst(src.rows, src.cols, src.type(), Scalar(0));
	Mat cvlt_mat(N, N, CV_64FC1, (double*)model_convolution);
	filter2D(src, dst, src.depth(), cvlt_mat);
	#endif

	#ifdef my_convolution
	int a = (N - 1) / 2;
	Mat tempImg(src.rows + 2*(N - 1), src.cols + 2 * (N - 1), src.type(), Scalar(0));
	Mat tempdst = (Mat_<double>(tempImg.rows - 2 * a, tempImg.cols - 2 * a));
	int ROI_x = (src.cols + 2 * (N - 1)) / 2 - src.cols / 2;//ROI矩形左上角的x坐标
	int ROI_y = (src.rows + 2 * (N - 1)) / 2 - src.rows / 2;//ROI矩形左上角的y坐标
	Rect ROIRect(ROI_x, ROI_y, src.cols, src.rows);//ROI矩形
	Mat tempImgROI2(tempImg, ROIRect);//tempImg的中间部分
	src.copyTo(tempImgROI2);//将原图复制到tempImg的中心
	for (int i = a; i < src.rows + 3 * a; i++)
	{
		for (int j = a; j < src.cols + 3 * a; j++)
		{
			sumLaplacian = 0;
			for (int k = -a; k <= a; k++)
			{
				for (int m = -a; m <= a; m++)
				{
					// 计算图像卷积
					sumLaplacian += tempImg.at<uchar>(i + k, j + m) *
					model_convolution[a + k][a + m];
				}
			}
			// 生成Lapace结果
			tempdst.at<double>(i-a, j-a) = sumLaplacian;
			//cout << tempdst.at<double>(i - a, j - a) << " ";//命令窗户输出
		}
		//cout << endl;//命令窗户输出
	}
	////截取原图大小
	ROI_x = tempdst.cols / 2 - src.cols / 2;//ROI矩形左上角的x坐标
	ROI_y = tempdst.rows / 2 - src.rows / 2;//ROI矩形左上角的y坐标
	Rect ROIRect2(ROI_x, ROI_y, src.cols, src.rows);//ROI矩形
	Mat dst(tempdst, ROIRect2);//tempdst的中间部分
	#endif

	//显示两种卷积的结果
	namedWindow(str, WINDOW_NORMAL);
	imshow(str, dst);

	//显示卷积前和卷积后的大小
	cout << src.rows << " * " << src.cols << endl;
	cout << dst.rows << " * " << dst.cols << endl;

	//二值化操作后的结果
	Mat result(src.rows, src.cols, src.type(), Scalar(0));

	//二值化操作
	//只设置阈值进行二值化
	#ifdef Only_threshold
	for (int i = 0; i < dst.rows - 1; i++)
	{
		for (int j = 0; j < dst.cols - 1; j++)
		{
			#ifdef filter2D_convolution
				result.at<uchar>(i, j) = dst.at<uchar>(i, j) > _threshold ? 255 : 0;
			#else
				result.at<uchar>(i, j) = dst.at<double>(i, j) > _threshold ? 255 : 0;
			#endif
		}
	}
	namedWindow(str + "二值化图像_阈值" + to_string(_threshold), WINDOW_NORMAL);
	imshow(str + "二值化图像_阈值" + to_string(_threshold), result);
	#endif
	//设置阈值并判断4邻域是否存在负值
	#ifdef threshold_neighborhood
	//3x3的邻域内，如果中心像素大于阈值，切周围存在负像素值，则此点为0交叉点
	for (int i = 0; i < dst.rows - 2; i++)
	{
		for (int j = 0; j < dst.cols - 2; j++)
		{
			result.at<uchar>(i + 1, j + 1) = 0;

			#ifdef filter2D_convolution
			if (dst.at<uchar>(i + 1, j + 1) >= _threshold && (
				dst.at<uchar>(i, j) <= 0 ||
				dst.at<uchar>(i, j + 1) <= 0 ||
				dst.at<uchar>(i + 1, j) <= 0 ||
				dst.at<uchar>(i + 2, j + 2) <= 0 ||
				dst.at<uchar>(i + 1, j + 2) <= 0 ||
				dst.at<uchar>(i + 2, j + 1) <= 0 ||
				dst.at<uchar>(i + 2, j) <= 0 ||
				dst.at<uchar>(i, j + 2) <= 0))
			{
				result.at<uchar>(i + 1, j + 1) = 255;
			}
			#else
			if (dst.at<double>(i + 1, j + 1) >= _threshold &&
				(dst.at<double>(i, j) <= -0.1 ||
				dst.at<double>(i, j + 1) <= -0.1 ||
				dst.at<double>(i + 1, j) <= -0.1 ||
				dst.at<double>(i + 2, j + 2) <= -0.1 ||
				dst.at<double>(i + 1, j + 2) <= -0.1 ||
				dst.at<double>(i + 2, j + 1) <= -0.1 ||
				dst.at<double>(i + 2, j) <= -0.1 ||
				dst.at<double>(i, j + 2) <= -0.1))
			{
				result.at<uchar>(i + 1, j + 1) = 255;
			}
			#endif

		}
	}
	namedWindow(str + "二值化_存在负邻域_阈值" + to_string(_threshold), WINDOW_NORMAL);
	imshow(str + "二值化_存在负邻域_阈值" + to_string(_threshold), result);
	#endif
	//判断4邻域乘积是否存在负值
	#ifdef neighborhood_negative_1
	//以P为中心的一个3*3领域，p点处的零交叉意味着至少有两个相对的领域像素的符号不同。
	//有四种要测试的情况：左/右、上/下，和两个对角。
	//P的值介于两个符号不同像素的值之间，就是P的值要小于两个符号不同的像素的绝对值。
	for (int y = 1; y < result.rows - 1; ++y)
	{
		for (int x = 1; x < result.cols - 1; ++x)
		{
			result.at<uchar>(y, x) = 0;
			#ifdef filter2D_convolution
			if (dst.at<uchar>(y - 1, x) *
				dst.at<uchar>(y + 1, x) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (dst.at<uchar>(y, x - 1) *
				dst.at<uchar>(y, x + 1) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (dst.at<uchar>(y + 1, x - 1) *
				dst.at<uchar>(y - 1, x + 1) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (dst.at<uchar>(y - 1, x - 1) *
				dst.at<uchar>(y + 1, x + 1) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			#else
			if (labs(dst.at<double>(y, x)) < labs(dst.at<double>(y - 1, x)) &&
				labs(dst.at<double>(y, x)) < labs(dst.at<double>(y + 1, x)) &&
				(dst.at<double>(y - 1, x) * dst.at<double>(y + 1, x)) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (labs(dst.at<double>(y, x)) < labs(dst.at<double>(y , x - 1)) &&
				labs(dst.at<double>(y, x)) < labs(dst.at<double>(y , x + 1)) &&
				(dst.at<double>(y , x- 1) * dst.at<double>(y , x+ 1)) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (labs(dst.at<double>(y, x)) < labs(dst.at<double>(y + 1, x - 1)) &&
				labs(dst.at<double>(y, x)) < labs(dst.at<double>(y - 1, x + 1)) &&
				(dst.at<double>(y + 1, x - 1) * dst.at<double>(y - 1, x + 1)) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (labs(dst.at<double>(y, x)) < labs(dst.at<double>(y - 1, x - 1)) &&
				labs(dst.at<double>(y, x)) < labs(dst.at<double>(y + 1, x + 1)) &&
				(dst.at<double>(y - 1, x - 1) * dst.at<double>(y + 1, x + 1)) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			#endif
		}
	}
	namedWindow(str + "二值化_4邻域乘积为负_1", WINDOW_NORMAL);
	imshow(str + "二值化_4邻域乘积为负_1", result);
	#endif
	//判断4邻域乘积是否存在负值
    #ifdef neighborhood_negative_2
	//以P为中心的一个3*3领域，p点处的零交叉意味着至少有两个相对的领域像素的符号不同。
	//有四种要测试的情况：左/右、上/下，和两个对角。
	//在这个基础上还要判断中心点p的绝对值是否是最接近0的点。
	for (int y = 1; y < result.rows - 1; ++y)
	{
		for (int x = 1; x < result.cols - 1; ++x)
		{
			result.at<uchar>(y, x) = 0;
			#ifdef filter2D_convolution
			if (dst.at<uchar>(y - 1, x) *
				dst.at<uchar>(y + 1, x) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (dst.at<uchar>(y, x - 1) *
				dst.at<uchar>(y, x + 1) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (dst.at<uchar>(y + 1, x - 1) *
				dst.at<uchar>(y - 1, x + 1) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			if (dst.at<uchar>(y - 1, x - 1) *
				dst.at<uchar>(y + 1, x + 1) <= 0)
			{
				result.at<uchar>(y, x) = 255;
			}
			#else
			//求4邻域内的最小值
			double min_dst_at = labs(dst.at<double>(y - 1, x - 1));
			for (int i = -1; i < 2; i++)
			{
				for (int j = -1; j < 2; j++)
				{
					if (labs(dst.at<double>(y + i, x + j)) <= min_dst_at)
					{
						min_dst_at = labs(dst.at<double>(y + i, x + j));
					}
				}
			}
			if (min_dst_at >= labs(dst.at<double>(y , x )))
			{
				if ((dst.at<double>(y - 1, x) * dst.at<double>(y + 1, x)) <= 0)
				{
					result.at<uchar>(y, x) = 255;
				}
				if ((dst.at<double>(y, x - 1) * dst.at<double>(y, x + 1)) <= 0)
				{
					result.at<uchar>(y, x) = 255;
				}
				if ((dst.at<double>(y + 1, x - 1) * dst.at<double>(y - 1, x + 1)) <= 0)
				{
					result.at<uchar>(y, x) = 255;
				}
				if ((dst.at<double>(y - 1, x - 1) * dst.at<double>(y + 1, x + 1)) <= 0)
				{
					result.at<uchar>(y, x) = 255;
				}

			}

			#endif
		}
	}
	namedWindow(str + "二值化_4邻域乘积为负_2", WINDOW_NORMAL);
	imshow(str + "二值化_4邻域乘积为负_2", result);
	#endif

	#ifdef store_result
	store_data(dst, dst_str_1);
	cout << "Store Done!" << endl;
	#endif

	return 0;
}
int main()
{
	Mat src = imread("Chapter3_1.pgm", 0);
	if (!src.data)
	{
		cout << "Error: read image" << endl;
		return -1;
	}
	cout << "src channel is: " << src.channels() << endl;
	namedWindow("src", WINDOW_NORMAL);
	imshow("src", src);

	#ifdef neighborhood_negative
	getLoGModel(dst_str_1, src, Sigma, 60);
	#else
	getLoGModel(dst_str_1, src, Sigma, 60);
	getLoGModel(dst_str_1, src, Sigma, 40);
	getLoGModel(dst_str_1, src, Sigma, 20);
	#endif

	waitKey(0);
	return 0;
}