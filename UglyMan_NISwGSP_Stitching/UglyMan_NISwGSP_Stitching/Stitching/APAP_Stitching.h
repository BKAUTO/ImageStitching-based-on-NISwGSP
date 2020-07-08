//
//  APAP_Stitching.h
//  UglyMan_Stitching
//
//  Created by uglyman.nothinglo on 2015/8/15.
//  Copyright (c) 2015 nothinglo. All rights reserved.
//

#ifndef __UglyMan_Stitiching__APAP_Stitching__
#define __UglyMan_Stitiching__APAP_Stitching__

#include "../Configure.h"
#include "../Util/Transform.h"
#include "../Feature/MultiImages.h"

class APAP_Stitching {
public:
    static void apap_project(const vector<Point2> & _p_src,
                             const vector<Point2> & _p_dst,
                             const vector<Point2> & _src,
                             vector<Point2>       & _dst,
                             vector<Mat>          & _homographies);

    static void apap_project(const vector<Point2> & _p_src,
                             const vector<Point2> & _p_dst,
                             const vector<KeyLine>& _l_src,
                             const vector<KeyLine>& _l_dst,
                             const vector<Point2> & _src,
                             vector<Point2>       & _dst,
                             vector<Mat>          & _homographies);

    static double distance(const Point2 P, const Point2 A, const Point2 B);
private:
};

#endif /* defined(__UglyMan_Stitiching__APAP_Stitching__) */
