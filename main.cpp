#include <iostream>
#include <string>
#include "src/imgproc/imgproc.h"
#include "src/utils/utility.h"
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "src/core/core.h"

using namespace std;
using namespace wcv;
using namespace cv;

void main()
{
	Mat src = imread("2.jpg",1);
	//src.convertTo(src, CV_16S);

	Image img;
	img.from_cvmat(src);

	Image ds;
	//wcv::threshold(img, ds, 50);
	//wcv::cvtColorGray(img, ds);
	//wcv::equalize(img, ds);
	//wcv::translationTrans(img, ds, -10, -10,Scalar4i::all(255));
	//wcv::mirrorTrans(img, ds, false);
	//wcv::tranpose(img, ds);
	//wcv::resize(img, ds, Size4i(src.cols * 4, src.rows * 4));
	//wcv::rotate(img, ds, 30);
	//wcv::copymakeBoarder(img, Size4i(7, 7), wcv::MIRROR, ds);

	float dat[] = {
		-1.,0.,1.,
		-2.,0.,2.,
		-1.,0.,1. };
	float dat1[] = {
		1,2,3,4,5,6,7,8,9,
		2,3,4,5,6,7,8,9,10,
		3,4,5,6,7,8,9,10,11,
		4,5,6,7,8,9,10,11,12,
		5,6,7,8,9,10,11,12,13,
		6,7,8,9,10,11,12,13,14,
		7,8,9,10,11,12,13,14,15,
		8,9,10,11,12,13,14,15,16,
		9,10,11,12,13,14,15,16,17
	};
	Mat32f k = Mat32f(3, 3, 3, (uchar*)(&dat1)),kk;
	//wcv::templateOp(img, k, ds, SAME);
	//Scalar4f kzz = k.dot(k);
	Mat32f lk = k.hadamardProduct(k);
	cout << lk.toString(false) << endl;

	vector<Mat8u> mvs;
	wcv::split(img, mvs);
	Mat8u imgg;
	wcv::merge(mvs, imgg);
	cv::Mat m0 = imgg.to_cvmat();


	//cv::Mat dst;
	//thresh(src, dst,20);
	Mat8u dss = img.subMat(Range4i(0, 100), Range4i(0, 100));
	Mat8u dsss;
	dsss= ds;

	Mat dst = dss.to_cvmat();
	imwrite("te.jpg", dst);
}

