//
//  ImageData.h
//  UglyMan_Stitching
//
//  Created by uglyman.nothinglo on 2015/8/15.
//  Copyright (c) 2015 nothinglo. All rights reserved.
//

#ifndef __UglyMan_Stitiching__ImageData__
#define __UglyMan_Stitiching__ImageData__

#include <memory>
#include "../Util/Statistics.h"
#include "../Feature/FeatureController.h"
#include "../Mesh/MeshGrid.h"

#include <iostream>
#include <opencv2/opencv_modules.hpp>
#include <opencv2/core.hpp>
#include <opencv2/line_descriptor.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/highgui.hpp>

#define MATCHES_DIST_THRESHOLD 25

using namespace cv::line_descriptor;

class LineData {
public:
    LineData(const Point2 & _a,
             const Point2 & _b,
             const double _width,
             const double _length);
    Point2 data[2];
    double width, length;
private:
};

typedef const bool (LINES_FILTER_FUNC)(const double _data, \
                                       const Statistics & _statistics);

LINES_FILTER_FUNC LINES_FILTER_NONE;
LINES_FILTER_FUNC LINES_FILTER_WIDTH;
LINES_FILTER_FUNC LINES_FILTER_LENGTH;


class ImageData {//单张图对象
public:
    string file_name, file_extension;
    const string * file_dir, * debug_dir;
    ImageData(const string & _file_dir,
              const string & _file_full_name,
              LINES_FILTER_FUNC * _width_filter,
              LINES_FILTER_FUNC * _length_filter,
              const string * _debug_dir = NULL);
    
    const Mat & getGreyImage() const;
    const vector<LineData> & getLines() const;
    const vector<Point2> & getFeaturePoints() const;
    const vector<FeatureDescriptor> & getFeatureDescriptors() const;
    
    const vector<KeyLine> & getKeyLines() const;
    const Mat & getLineFeatureDescriptor() const;

    void clear();
    
    Mat img, rgba_img, alpha_mask;
    unique_ptr<Mesh2D> mesh_2d;
    
    mutable vector<KeyLine> keylines;

private:
    LINES_FILTER_FUNC * width_filter, * length_filter;
    
    mutable Mat grey_img;
    mutable vector<LineData> img_lines;
    mutable vector<Point2> feature_points;
    mutable vector<FeatureDescriptor> feature_descriptors;

    //线特征匹配所需参数
    mutable Mat descr;
};

#endif /* defined(__UglyMan_Stitiching__ImageData__) */
