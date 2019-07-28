//opencv版本:OpenCV2.7.13
//VS版本:VS2015
//Author:zqs
//time:2018-0322
//name:turn picture

#include<iostream>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<math.h>
#include<fstream>
#include<string>

using namespace std;
using namespace cv;

#define SCALE 1   //缩放比例
#define pi 3.1415926

/*
//interpolation seting,if set '0' it is nearest neighbor interpolation,if set '1' it is bilinear interpolation,
//if set '2' it is bicubic interpolation,if set '3' it is area-based (or super) interpolation,
//if set '4' it is Lanczos interpolation over 8x8 neighborhood.
*/
#define interpolation 1

Mat turn_picture(Mat src,int angle,int _interpolation)
{
	Mat dst;//输出旋转后的图像
	double angle2 = angle * pi / 180;
	double w = 0., h = 0., w_r = 0., h_r = 0.;

	h = src.rows;
	w = src.cols;
	w_r = w*cos(angle2) + h*sin(angle2);//输出图像的宽度或高度
	h_r = h*cos(angle2) + w*sin(angle2);//为了保证无论如何旋转都能放下

	//建立临时图像，长宽都是源图像的对角线长度，将源图像复制到临时图像的中心后再变换
	Mat tempImg(w_r, h_r, src.type(), Scalar(0));
	//临时图像，大小和输出图像一样大
	int ROI_x = w_r / 2 - src.cols / 2;//ROI矩形左上角的x坐标
	int ROI_y = h_r / 2 - src.rows / 2;//ROI矩形左上角的y坐标
	Rect ROIRect(ROI_x, ROI_y, src.cols, src.rows);//ROI矩形
	Mat tempImgROI2(tempImg, ROIRect);//tempImg的中间部分
	src.copyTo(tempImgROI2);//将原图复制到tempImg的中心

	Point2f center(w_r / 2, h_r / 2);//旋转中心
	Mat M = getRotationMatrix2D(center, angle, SCALE);//计算旋转的仿射变换矩阵

													  //输出看看算出的矩阵是什么
	cout << "变换矩阵：" << endl;
	cout << M.at<double>(0, 0) << "," << M.at<double>(0, 1) << "," << M.at<double>(0, 2) << "," << endl;
	cout << M.at<double>(1, 0) << "," << M.at<double>(1, 1) << "," << M.at<double>(1, 2) << "," << endl;

	warpAffine(tempImg, dst, M, Size(w_r, h_r), interpolation);//仿射变换

	return dst;
}
Mat turn_picture_2(Mat src_src ,Mat src, int angle, int _interpolation)
{
	Mat dst_2;//输出旋转后的图像
	double w = 0., h = 0.;

	h = src.rows;
	w = src.cols;

	Point2f center(w / 2, h / 2);//旋转中心
	Mat M = getRotationMatrix2D(center, angle, SCALE);//计算旋转的仿射变换矩阵

													  //输出看看算出的矩阵是什么
	cout << "变换矩阵：" << endl;
	cout << M.at<double>(0, 0) << "," << M.at<double>(0, 1) << "," << M.at<double>(0, 2) << "," << endl;
	cout << M.at<double>(1, 0) << "," << M.at<double>(1, 1) << "," << M.at<double>(1, 2) << "," << endl;

	warpAffine(src, dst_2, M, Size(w, h), interpolation);//仿射变换

	int ROI_x = w / 2 - src_src.cols / 2;//ROI矩形左上角的x坐标
	int ROI_y = h / 2 - src_src.rows / 2;//ROI矩形左上角的y坐标
	Rect ROIRect(ROI_x, ROI_y, src_src.cols, src_src.rows);//ROI矩形
	Mat tempImgROI2(dst_2, ROIRect);//dst_2的中间部分
	cout << tempImgROI2.cols << endl;
	cout << tempImgROI2.rows << endl;
	return tempImgROI2;
}

void store_data(Mat &image, string txt_name)
{
	ofstream outfile;
	outfile.open(txt_name + ".txt"); //存放数据的文件名

	Mat input_image = image.clone();//复制实参到临时变量
	int n_r = input_image.rows; // 行数
	//列数*通道数等于每一行元素的个数
	int n_c = input_image.cols * input_image.channels();
	for (int j = 0; j<n_r; j++)
	{
		uchar* data = input_image.ptr<uchar>(j);
		for (int i = 0; i<n_c; i++)
		{
			if (outfile.is_open())
			{
				outfile << (int)data[i] << ","; //程序中处理的数据
				//outfile.close();
			}
			else
			{
				cout << "不能打开文件!" << endl;
			}
			//cout << (int)data[i] << ",";
		}
		if (outfile.is_open())
		{
			outfile << endl; //message是程序中处理的数据
			//outfile.close();
		}
		else
		{
			cout << "不能打开文件!" << endl;
		}
		//cout << endl;
	}
	outfile.close();
}

int main()
{
	Mat src = imread("Chapter2_1.pgm", 0);
	string txt_src = "Chapter2_1_interpolation_" + to_string(interpolation);
	string txt_turn_picture = "turn_picture_interpolation_" + to_string(interpolation);

	cout << src.cols << endl;
	cout << src.rows << endl;

	store_data(src, txt_src);

	int angle_1 = 15;//旋转角度(正值表示逆时针旋转)
	int angle_2 = -15;
	//显示
	imshow("src", src);
	imshow("turn picture", turn_picture(src,angle_1,interpolation));

	imshow("turn picture_2", turn_picture_2(src,turn_picture(src, angle_1, interpolation), angle_2, interpolation));

	store_data(turn_picture_2(src, turn_picture(src, angle_1, interpolation), angle_2, interpolation), txt_turn_picture);

	waitKey(0);
	return 0;
}