//
//  APAP_Stitching.cpp
//  UglyMan_Stitching
//
//  Created by uglyman.nothinglo on 2015/8/15.
//  Copyright (c) 2015 nothinglo. All rights reserved.
//

#include "APAP_Stitching.h"

double APAP_Stitching::distance(const Point2 P, const Point2 A, const Point2 B){
    //计算r |AB| |AP| |BP| |PC|
	
	double ab =	sqrt(pow((B.x-A.x),2) +pow((B.y-A.y),2)); // |AB|
	double ap = sqrt(pow((P.x-A.x),2) +pow((P.y-A.y),2)); // |AP|
	double bp = sqrt(pow((P.x-B.x),2) +pow((P.y-B.y),2)); // |BP|
	double r = 0;
	if (ab > 0)
	{
		r = ((P.x - A.x)*(B.x - A.x) + (P.y - A.y)*(B.y - A.y)) / pow(ab, 2);
	} //r 
	else
	{
		cout << "no lines" << endl;
	}

	//double distance = 0;
	double distance = 0;
	if(ab>0)
	{
		if(r>=1)
			distance = bp;
		else if (r>0&&r<1)
			distance = sqrt(pow(ap,2)-r*r*pow(ab,2));
		else
			distance = ap;
	}

	//显示结果,可以带入点进行测试
    /*
	cout<<"|r|"<<r<<endl;	
	cout<<"|AP|"<<ap<<endl;
	cout<<"|BP|"<<bp<<endl;
	cout<<"|AB|"<<ab<<endl;
	cout<<"distance is : "<<distance<<endl;
	*/
	return distance;
}


void APAP_Stitching::apap_project(const vector<Point2> & _p_src,
                                  const vector<Point2> & _p_dst,
                                  const vector<Point2> & _src,
                                  vector<Point2>       & _dst,
                                  vector<Mat>          & _homographies) {
    vector<Point2> nf1, nf2, cf1, cf2;
    Mat N1, N2, C1, C2;
    N1 = getNormalize2DPts(_p_src, nf1);
    N2 = getNormalize2DPts(_p_dst, nf2);
    C1 = getConditionerFromPts(nf1);
    C2 = getConditionerFromPts(nf2);
    cf1.reserve(nf1.size());
    cf2.reserve(nf2.size());
    for(int i = 0; i < nf1.size(); ++i) {
        cf1.emplace_back(nf1[i].x * C1.at<double>(0, 0) + C1.at<double>(0, 2),
                         nf1[i].y * C1.at<double>(1, 1) + C1.at<double>(1, 2));

        cf2.emplace_back(nf2[i].x * C2.at<double>(0, 0) + C2.at<double>(0, 2),
                         nf2[i].y * C2.at<double>(1, 1) + C2.at<double>(1, 2));
    }
    double sigma_inv_2 = 1. / (APAP_SIGMA * APAP_SIGMA), gamma = APAP_GAMMA;
    MatrixXd A = MatrixXd::Zero(cf1.size() * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);//矩阵A的大小为2N×9
    
#ifndef NDEBUG
    if(_dst.empty() == false) {
        _dst.clear();
        printError("F(apap_project) dst is not empty");
    }
    if(_homographies.empty() == false) {
        _homographies.clear();
        printError("F(apap_project) homographies is not empty");
    }
#endif
    _dst.reserve(_src.size());
    _homographies.reserve(_src.size());
    for(int i = 0; i < _src.size(); ++i) {
        for(int j = 0; j < _p_src.size(); ++j) {
            Point2 d = _src[i] - _p_src[j];
            double www = MAX(gamma, exp(-sqrt(d.x * d.x + d.y * d.y) * sigma_inv_2));
            A(2*j  , 0) = www * cf1[j].x;
            A(2*j  , 1) = www * cf1[j].y;
            A(2*j  , 2) = www * 1;
            A(2*j  , 6) = www * -cf2[j].x * cf1[j].x;
            A(2*j  , 7) = www * -cf2[j].x * cf1[j].y;
            A(2*j  , 8) = www * -cf2[j].x;
            
            A(2*j+1, 3) = www * cf1[j].x;
            A(2*j+1, 4) = www * cf1[j].y;
            A(2*j+1, 5) = www * 1;
            A(2*j+1, 6) = www * -cf2[j].y * cf1[j].x;
            A(2*j+1, 7) = www * -cf2[j].y * cf1[j].y;
            A(2*j+1, 8) = www * -cf2[j].y;
        }
        JacobiSVD<MatrixXd, HouseholderQRPreconditioner> jacobi_svd(A, ComputeThinV);
        MatrixXd V = jacobi_svd.matrixV();
        Mat H(3, 3, CV_64FC1);
        for(int j = 0; j < V.rows(); ++j) {
            H.at<double>(j / 3, j % 3) = V(j, V.rows() - 1);
        }
        H = C2.inv() * H * C1;
        H = N2.inv() * H * N1;

        //std::cout << "local homography H " << i << " = "<< H << std::endl;

        _dst.emplace_back(applyTransform3x3(_src[i].x, _src[i].y, H));
        _homographies.emplace_back(H);
    }
}

void APAP_Stitching::apap_project(const vector<Point2> & _p_src,
                                  const vector<Point2> & _p_dst,
                                  const vector<KeyLine>& _l_src,
                                  const vector<KeyLine>& _l_dst,
                                  const vector<Point2> & _src,
                                  vector<Point2>       & _dst,
                                  vector<Mat>          & _homographies) {
    vector<Point2> lineStart_src, lineEnd_src, lineMid_src, lineQuar1_src, lineQuar2_src, lineStart_dst, lineEnd_dst, lineMid_dst, lineQuar1_dst, lineQuar2_dst;
    
    lineStart_src.reserve(_l_src.size());
    lineEnd_src.reserve(_l_src.size());
    lineMid_src.reserve(_l_src.size());
    lineStart_dst.reserve(_l_src.size());
    lineEnd_dst.reserve(_l_src.size());
    lineMid_dst.reserve(_l_src.size());
    lineQuar1_src.reserve(_l_src.size());
    lineQuar2_src.reserve(_l_src.size());
    lineQuar1_dst.reserve(_l_src.size());
    lineQuar2_dst.reserve(_l_src.size());

    for(int i = 0; i < _l_src.size(); i++){
        lineStart_src.emplace_back(_l_src[i].getStartPoint()); 
        lineEnd_src.emplace_back(_l_src[i].getEndPoint());
        lineMid_src.emplace_back((_l_src[i].getStartPoint().x + _l_src[i].getEndPoint().x) / 2.0, 
                                 (_l_src[i].getStartPoint().y + _l_src[i].getEndPoint().y) / 2.0);
        lineQuar1_src.emplace_back((_l_src[i].getStartPoint().x + lineMid_src[i].x) / 2.0, 
                                   (_l_src[i].getStartPoint().y + lineMid_src[i].y) / 2.0);
        lineQuar2_src.emplace_back((_l_src[i].getEndPoint().x + lineMid_src[i].x) / 2.0, 
                                   (_l_src[i].getEndPoint().y + lineMid_src[i].y) / 2.0);

        lineStart_dst.emplace_back(_l_dst[i].getStartPoint()); 
        lineEnd_dst.emplace_back(_l_dst[i].getEndPoint());
        lineMid_dst.emplace_back((_l_dst[i].getStartPoint().x + _l_dst[i].getEndPoint().x) / 2.0, 
                                 (_l_dst[i].getStartPoint().y + _l_dst[i].getEndPoint().y) / 2.0);
        lineQuar1_dst.emplace_back((_l_dst[i].getStartPoint().x + lineMid_dst[i].x) / 2.0, 
                                   (_l_dst[i].getStartPoint().y + lineMid_dst[i].y) / 2.0);
        lineQuar2_dst.emplace_back((_l_dst[i].getEndPoint().x + lineMid_dst[i].x) / 2.0, 
                                   (_l_dst[i].getEndPoint().y + lineMid_dst[i].y) / 2.0);
    }
    /*
    std::cout << "point_size : " << _p_src.size() << std::endl;
    std::cout << "line_size : " << _l_src.size() << std::endl;
    std::cout << "lineStart_size : " << lineStart_src.size() << std::endl;
    std::cout << "lineEnd_size : " << lineEnd_src.size() << std::endl;
    std::cout << "lineMid_size : " << lineMid_src.size() << std::endl;
    std::cout << "lineQuar1_size : " << lineQuar1_dst.size() << std::endl;
    std::cout << "lineQuar2_size : " << lineQuar2_dst.size() << std::endl;
    */
    vector<Point2> nf1, nf2, cf1, cf2;
    Mat N1, N2, C1, C2;

    vector<Point2> allPoint_src, allPoint_dst;
    allPoint_src.reserve(_p_src.size() + 5 * _l_src.size());
    allPoint_dst.reserve(_p_src.size() + 5 * _l_src.size());
    allPoint_src.insert(allPoint_src.end(), _p_src.begin(), _p_src.end());
    allPoint_src.insert(allPoint_src.end(), lineStart_src.begin(), lineStart_src.end());
    allPoint_src.insert(allPoint_src.end(), lineEnd_src.begin(), lineEnd_src.end());
    allPoint_src.insert(allPoint_src.end(), lineMid_src.begin(), lineMid_src.end());
    allPoint_src.insert(allPoint_src.end(), lineQuar1_src.begin(), lineQuar1_src.end());
    allPoint_src.insert(allPoint_src.end(), lineQuar2_src.begin(), lineQuar2_src.end());

    allPoint_dst.insert(allPoint_dst.end(), _p_dst.begin(), _p_dst.end());
    allPoint_dst.insert(allPoint_dst.end(), lineStart_dst.begin(), lineStart_dst.end());
    allPoint_dst.insert(allPoint_dst.end(), lineEnd_dst.begin(), lineEnd_dst.end());
    allPoint_dst.insert(allPoint_dst.end(), lineMid_dst.begin(), lineMid_dst.end());
    allPoint_dst.insert(allPoint_dst.end(), lineQuar1_dst.begin(), lineQuar1_dst.end());
    allPoint_dst.insert(allPoint_dst.end(), lineQuar2_dst.begin(), lineQuar2_dst.end());

    N1 = getNormalize2DPts(allPoint_src, nf1);
    N2 = getNormalize2DPts(allPoint_dst, nf2);
    C1 = getConditionerFromPts(nf1);
    C2 = getConditionerFromPts(nf2);
    cf1.reserve(nf1.size());
    cf2.reserve(nf2.size());
    for(int i = 0; i < nf1.size(); ++i) {
        cf1.emplace_back(nf1[i].x * C1.at<double>(0, 0) + C1.at<double>(0, 2),
                         nf1[i].y * C1.at<double>(1, 1) + C1.at<double>(1, 2));

        cf2.emplace_back(nf2[i].x * C2.at<double>(0, 0) + C2.at<double>(0, 2),
                         nf2[i].y * C2.at<double>(1, 1) + C2.at<double>(1, 2));
    }

    double sigma_inv_2 = 1. / (APAP_SIGMA * APAP_SIGMA), gamma = APAP_GAMMA;
    MatrixXd A = MatrixXd::Zero(_p_src.size() * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);//矩阵A的大小为2N×9
    MatrixXd B1 = MatrixXd::Zero(_l_src.size() * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);//矩阵B的大小为2M×9
    MatrixXd B2 = MatrixXd::Zero(_l_src.size() * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);//矩阵B的大小为2M×9
    MatrixXd B3 = MatrixXd::Zero(_l_src.size() * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);//矩阵B的大小为2M×9
    MatrixXd B4 = MatrixXd::Zero(_l_src.size() * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);//矩阵B的大小为2M×9
    MatrixXd B5 = MatrixXd::Zero(_l_src.size() * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);//矩阵B的大小为2M×9
    MatrixXd C = MatrixXd::Zero((_p_src.size() + 5 * _l_src.size()) * DIMENSION_2D,
                                HOMOGRAPHY_VARIABLES_COUNT);

#ifndef NDEBUG
    if(_dst.empty() == false) {
        _dst.clear();
        printError("F(apap_project) dst is not empty");
    }
    if(_homographies.empty() == false) {
        _homographies.clear();
        printError("F(apap_project) homographies is not empty");
    }
#endif
    _dst.reserve(_src.size());
    _homographies.reserve(_src.size());
    for(int i = 0; i < _src.size(); ++i) {
        for(int j = 0; j < _p_src.size(); ++j) {
            Point2 d = _src[i] - _p_src[j];
            double www = MAX(gamma, exp(-sqrt(d.x * d.x + d.y * d.y) * sigma_inv_2));
            //std::cout << "www : " << www << std::endl;
            A(2*j  , 0) = www * cf1[j].x;
            A(2*j  , 1) = www * cf1[j].y;
            A(2*j  , 2) = www * 1;
            A(2*j  , 6) = www * -cf2[j].x * cf1[j].x;
            A(2*j  , 7) = www * -cf2[j].x * cf1[j].y;
            A(2*j  , 8) = www * -cf2[j].x;
            
            A(2*j+1, 3) = www * cf1[j].x;
            A(2*j+1, 4) = www * cf1[j].y;
            A(2*j+1, 5) = www * 1;
            A(2*j+1, 6) = www * -cf2[j].y * cf1[j].x;
            A(2*j+1, 7) = www * -cf2[j].y * cf1[j].y;
            A(2*j+1, 8) = www * -cf2[j].y;
        }
        
        for(int k = 0; k < _l_src.size(); ++k){
            int line_idx1 = k + _p_src.size();
            int line_idx2 = line_idx1 + _l_src.size();
            int line_idx3 = line_idx2 + _l_src.size();
            int line_idx4 = line_idx3 + _l_src.size();
            int line_idx5 = line_idx4 + _l_src.size();
            double wwl = MAX(gamma, exp(-distance(_src[i], _l_src[k].getStartPoint(), _l_src[k].getEndPoint()) * sigma_inv_2));
            //std::cout << "wwl : " << wwl << std::endl;
            B1(2*k  , 0) = wwl * cf1[line_idx1].x;
            B1(2*k  , 1) = wwl * cf1[line_idx1].y;
            B1(2*k  , 2) = wwl * 1;
            B1(2*k  , 6) = wwl * -cf2[line_idx1].x * cf1[line_idx1].x;
            B1(2*k  , 7) = wwl * -cf2[line_idx1].x * cf1[line_idx1].y;
            B1(2*k  , 8) = wwl * -cf2[line_idx1].x;
            
            B1(2*k+1, 3) = wwl * cf1[line_idx1].x;
            B1(2*k+1, 4) = wwl * cf1[line_idx1].y;
            B1(2*k+1, 5) = wwl * 1;
            B1(2*k+1, 6) = wwl * -cf2[line_idx1].y * cf1[line_idx1].x;
            B1(2*k+1, 7) = wwl * -cf2[line_idx1].y * cf1[line_idx1].y;
            B1(2*k+1, 8) = wwl * -cf2[line_idx1].y;

            B2(2*k  , 0) = wwl * cf1[line_idx2].x;
            B2(2*k  , 1) = wwl * cf1[line_idx2].y;
            B2(2*k  , 2) = wwl * 1;
            B2(2*k  , 6) = wwl * -cf2[line_idx2].x * cf1[line_idx2].x;
            B2(2*k  , 7) = wwl * -cf2[line_idx2].x * cf1[line_idx2].y;
            B2(2*k  , 8) = wwl * -cf2[line_idx2].x;
            
            B2(2*k+1, 3) = wwl * cf1[line_idx2].x;
            B2(2*k+1, 4) = wwl * cf1[line_idx2].y;
            B2(2*k+1, 5) = wwl * 1;
            B2(2*k+1, 6) = wwl * -cf2[line_idx2].y * cf1[line_idx2].x;
            B2(2*k+1, 7) = wwl * -cf2[line_idx2].y * cf1[line_idx2].y;
            B2(2*k+1, 8) = wwl * -cf2[line_idx2].y;

            B3(2*k  , 0) = wwl * cf1[line_idx3].x;
            B3(2*k  , 1) = wwl * cf1[line_idx3].y;
            B3(2*k  , 2) = wwl * 1;
            B3(2*k  , 6) = wwl * -cf2[line_idx3].x * cf1[line_idx3].x;
            B3(2*k  , 7) = wwl * -cf2[line_idx3].x * cf1[line_idx3].y;
            B3(2*k  , 8) = wwl * -cf2[line_idx3].x;
            
            B3(2*k+1, 3) = wwl * cf1[line_idx3].x;
            B3(2*k+1, 4) = wwl * cf1[line_idx3].y;
            B3(2*k+1, 5) = wwl * 1;
            B3(2*k+1, 6) = wwl * -cf2[line_idx3].y * cf1[line_idx3].x;
            B3(2*k+1, 7) = wwl * -cf2[line_idx3].y * cf1[line_idx3].y;
            B3(2*k+1, 8) = wwl * -cf2[line_idx3].y;

            B4(2*k  , 0) = wwl * cf1[line_idx4].x;
            B4(2*k  , 1) = wwl * cf1[line_idx4].y;
            B4(2*k  , 2) = wwl * 1;
            B4(2*k  , 6) = wwl * -cf2[line_idx4].x * cf1[line_idx4].x;
            B4(2*k  , 7) = wwl * -cf2[line_idx4].x * cf1[line_idx4].y;
            B4(2*k  , 8) = wwl * -cf2[line_idx4].x;
            
            B4(2*k+1, 3) = wwl * cf1[line_idx4].x;
            B4(2*k+1, 4) = wwl * cf1[line_idx4].y;
            B4(2*k+1, 5) = wwl * 1;
            B4(2*k+1, 6) = wwl * -cf2[line_idx4].y * cf1[line_idx4].x;
            B4(2*k+1, 7) = wwl * -cf2[line_idx4].y * cf1[line_idx4].y;
            B4(2*k+1, 8) = wwl * -cf2[line_idx4].y;

            B5(2*k  , 0) = wwl * cf1[line_idx5].x;
            B5(2*k  , 1) = wwl * cf1[line_idx5].y;
            B5(2*k  , 2) = wwl * 1;
            B5(2*k  , 6) = wwl * -cf2[line_idx5].x * cf1[line_idx5].x;
            B5(2*k  , 7) = wwl * -cf2[line_idx5].x * cf1[line_idx5].y;
            B5(2*k  , 8) = wwl * -cf2[line_idx5].x;
            
            B5(2*k+1, 3) = wwl * cf1[line_idx5].x;
            B5(2*k+1, 4) = wwl * cf1[line_idx5].y;
            B5(2*k+1, 5) = wwl * 1;
            B5(2*k+1, 6) = wwl * -cf2[line_idx5].y * cf1[line_idx5].x;
            B5(2*k+1, 7) = wwl * -cf2[line_idx5].y * cf1[line_idx5].y;
            B5(2*k+1, 8) = wwl * -cf2[line_idx5].y;
        }
        
        C << A,
             B1,
             B2,
             B3,
             B4,
             B5;
        
        JacobiSVD<MatrixXd, HouseholderQRPreconditioner> jacobi_svd(C, ComputeThinV);
        MatrixXd V = jacobi_svd.matrixV();
        Mat H(3, 3, CV_64FC1);
        for(int j = 0; j < V.rows(); ++j) {
            H.at<double>(j / 3, j % 3) = V(j, V.rows() - 1);
        }
        H = C2.inv() * H * C1;
        H = N2.inv() * H * N1;

        //std::cout << "local homography H " << i << " = "<< H << std::endl;

        _dst.emplace_back(applyTransform3x3(_src[i].x, _src[i].y, H));
        _homographies.emplace_back(H);
    }
}