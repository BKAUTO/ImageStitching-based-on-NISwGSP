//
//  MultiImages.cpp
//  UglyMan_Stitching
//
//  Created by uglyman.nothinglo on 2015/8/15.
//  Copyright (c) 2015 nothinglo. All rights reserved.
//

#include "MultiImages.h"
#include "opencv2/stitching.hpp"
#include "opencv2/stitching/detail/blenders.hpp"
using namespace cv::detail;
#include "opencv2/stitching/detail/util.hpp"
#include "opencv2/stitching/detail/exposure_compensate.hpp"

MultiImages::MultiImages(const string & _file_name,
                         LINES_FILTER_FUNC * _width_filter,
                         LINES_FILTER_FUNC * _length_filter) : parameter(_file_name) {
    
    for(int i = 0; i < parameter.image_file_full_names.size(); ++i) {
#ifndef NDEBUG
        images_data.emplace_back(parameter.file_dir,
                                 parameter.image_file_full_names[i],
                                 _width_filter,
                                 _length_filter,
                                 &parameter.debug_dir);
#else
        images_data.emplace_back(parameter.file_dir,
                                 parameter.image_file_full_names[i],
                                 _width_filter,
                                 _length_filter);
#endif
    }
}

void MultiImages::doFeatureMatching() const {
    const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
    
    images_features.resize(images_data.size());
    images_features_mask.resize(images_data.size());
    for(int i = 0; i < images_data.size(); ++i) {
        const vector<Point2> & vertices = images_data[i].mesh_2d->getVertices();
        images_features_mask[i].resize(vertices.size(), false);
        for(int j = 0; j < vertices.size(); ++j) {
            images_features[i].keypoints.emplace_back(vertices[j], 0);
        }
    }
    pairwise_matches.resize(images_data.size() * images_data.size());
    
    apap_homographies.resize(images_data.size());
    apap_overlap_mask.resize(images_data.size());
    apap_matching_points.resize(images_data.size());
    for(int i = 0; i < images_data.size(); ++i) {
        apap_homographies[i].resize(images_data.size());
        apap_overlap_mask[i].resize(images_data.size());
        apap_matching_points[i].resize(images_data.size());
    }
    const vector<vector<vector<Point2> > > & feature_matches = getFeatureMatches();
    const vector<vector<vector<KeyLine> > > & line_feature_matches = getLineFeatureMatches();
    
    for(int i = 0; i < images_match_graph_pair_list.size(); ++i) {
        const pair<int, int> & match_pair = images_match_graph_pair_list[i];
        const int & m1 = match_pair.first, & m2 = match_pair.second;

        vector<Point2> f1, f2;
        f1.reserve(line_feature_matches[m1][m2].size());
        f2.reserve(line_feature_matches[m1][m2].size());
        for(int i = 0; i < line_feature_matches[m1][m2].size(); i++){
            f1.emplace_back(line_feature_matches[m1][m2][i].getStartPoint());
            f2.emplace_back(line_feature_matches[m2][m1][i].getStartPoint());
        }
        Mat image_of_feauture_pairs = getImageOfFeaturePairs(images_data[match_pair.first ].img,
                                                         images_data[match_pair.second].img,
                                                         f1, f2);
        string _name = "TEST";
        imwrite(parameter.debug_dir +
            "feature_pairs-" + _name + "-" +
            images_data[match_pair.first ].file_name  + "-" +
            images_data[match_pair.second].file_name  + "-" +
            to_string(f1.size()) +
            images_data[match_pair.first ].file_extension, image_of_feauture_pairs);

        APAP_Stitching::apap_project(feature_matches[m1][m2],
                                     feature_matches[m2][m1],
                                     line_feature_matches[m1][m2],
                                     line_feature_matches[m2][m1],
                                     images_data[m1].mesh_2d->getVertices(), apap_matching_points[m1][m2], apap_homographies[m1][m2]);
        APAP_Stitching::apap_project(feature_matches[m2][m1],
                                     feature_matches[m1][m2],
                                     line_feature_matches[m2][m1],
                                     line_feature_matches[m1][m2],
                                     images_data[m2].mesh_2d->getVertices(), apap_matching_points[m2][m1], apap_homographies[m2][m1]);
        
        const int PAIR_SIZE = 2;
        const vector<Point2> * out_dst[PAIR_SIZE] = { &apap_matching_points[m1][m2], &apap_matching_points[m2][m1] };
        
        apap_overlap_mask[m1][m2].resize(apap_homographies[m1][m2].size(), false);
        apap_overlap_mask[m2][m1].resize(apap_homographies[m2][m1].size(), false);
        
        const int pm_index = m1 * (int)images_data.size() + m2;
        const int m_index[PAIR_SIZE] = {m2, m1};
        vector<DMatch> & D_matches = pairwise_matches[pm_index].matches;
        for(int j = 0; j < PAIR_SIZE; ++j) {
            for(int k = 0; k < out_dst[j]->size(); ++k) {
                if((*out_dst[j])[k].x >= 0 && (*out_dst[j])[k].y >= 0 &&
                   (*out_dst[j])[k].x <= images_data[m_index[j]].img.cols &&
                   (*out_dst[j])[k].y <= images_data[m_index[j]].img.rows) {
                    if(j) {
                        apap_overlap_mask[m2][m1][k] = true;
                        D_matches.emplace_back(images_features[m_index[j]].keypoints.size(), k, 0);
                        images_features_mask[m2][k] = true;
                    } else {
                        apap_overlap_mask[m1][m2][k] = true;
                        D_matches.emplace_back(k, images_features[m_index[j]].keypoints.size(), 0);
                        images_features_mask[m1][k] = true;
                    }
                    images_features[m_index[j]].keypoints.emplace_back((*out_dst[j])[k], 0);
                }
            }
        }
        pairwise_matches[pm_index].confidence  = 2.; /*** need > 1.f ***/
        pairwise_matches[pm_index].src_img_idx = m1;
        pairwise_matches[pm_index].dst_img_idx = m2;
        pairwise_matches[pm_index].inliers_mask.resize(D_matches.size(), 1);
        pairwise_matches[pm_index].num_inliers = (int)D_matches.size();
        pairwise_matches[pm_index].H = apap_homographies[m1][m2].front(); /*** for OpenCV findMaxSpanningTree funtion ***/
        
        Mat matchImage1;
        drawMatches(images_data[m1].img, images_features[m1].keypoints, images_data[m2].img,
			images_features[m2].keypoints, pairwise_matches[pm_index].matches, matchImage1);
		//namedWindow("match", WINDOW_NORMAL);
		//imshow("match", matchImage1);
		//waitKey();
        imwrite(parameter.debug_dir +
            "align_result-" + _name + "-" +
            images_data[match_pair.first ].file_name  + "-" +
            images_data[match_pair.second].file_name  + "-" +
            images_data[match_pair.first ].file_extension, matchImage1);
    }
}

const vector<detail::ImageFeatures> & MultiImages::getImagesFeaturesByMatchingPoints() const {
    if(images_features.empty()) {
        doFeatureMatching();
    }
    return images_features;
}

const vector<detail::MatchesInfo> & MultiImages::getPairwiseMatchesByMatchingPoints() const {
    if(pairwise_matches.empty()) {
        doFeatureMatching();
    }
    return pairwise_matches;
}

const vector<detail::CameraParams> & MultiImages::getCameraParams() const {
    if(camera_params.empty()) {
        camera_params.resize(images_data.size());
        /*** Focal Length ***/
        const vector<vector<vector<bool> > > & apap_overlap_mask = getAPAPOverlapMask();
        const vector<vector<vector<Mat> > >  & apap_homographies = getAPAPHomographies();
        
        vector<Mat> translation_matrix;
        translation_matrix.reserve(images_data.size());
        for(int i = 0; i < images_data.size(); ++i) {
            Mat T(3, 3, CV_64FC1);
            T.at<double>(0, 0) = T.at<double>(1, 1) = T.at<double>(2, 2) = 1;
            T.at<double>(0, 2) = images_data[i].img.cols * 0.5;
            T.at<double>(1, 2) = images_data[i].img.rows * 0.5;
            T.at<double>(0, 1) = T.at<double>(1, 0) = T.at<double>(2, 0) = T.at<double>(2, 1) = 0;
            translation_matrix.emplace_back(T);
        }
        vector<vector<double> > image_focal_candidates;
        image_focal_candidates.resize(images_data.size());
        for(int i = 0; i < images_data.size(); ++i) {
            for(int j = 0; j < images_data.size(); ++j) {
                for(int k = 0; k < apap_overlap_mask[i][j].size(); ++k) {
                    if(apap_overlap_mask[i][j][k]) {
                        double f0, f1;
                        bool f0_ok, f1_ok;
                        Mat H = translation_matrix[j].inv() * apap_homographies[i][j][k] * translation_matrix[i];
                        detail::focalsFromHomography(H / H.at<double>(2, 2),
                                                     f0, f1, f0_ok, f1_ok);
                        if(f0_ok && f1_ok) {
                            image_focal_candidates[i].emplace_back(f0);
                            image_focal_candidates[j].emplace_back(f1);
                        }
                    }
                }
            }
        }
        for(int i = 0; i < camera_params.size(); ++i) {
            if(image_focal_candidates[i].empty()) {
                camera_params[i].focal = images_data[i].img.cols + images_data[i].img.rows;
            } else {
                Statistics::getMedianWithoutCopyData(image_focal_candidates[i], camera_params[i].focal);
            }
        }
        /********************/
        /*** 3D Rotations ***/
        vector<vector<Mat> > relative_3D_rotations;
        relative_3D_rotations.resize(images_data.size());
        for(int i = 0; i < relative_3D_rotations.size(); ++i) {
            relative_3D_rotations[i].resize(images_data.size());
        }
        const vector<detail::ImageFeatures> & images_features    = getImagesFeaturesByMatchingPoints();
        const vector<detail::MatchesInfo>   & pairwise_matches   = getPairwiseMatchesByMatchingPoints();
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        for(int i = 0; i < images_match_graph_pair_list.size(); ++i) {
            const pair<int, int> & match_pair = images_match_graph_pair_list[i];
            const int & m1 = match_pair.first, & m2 = match_pair.second;
            const int m_index = m1 * (int)images_data.size() + m2;
            const detail::MatchesInfo & matches_info = pairwise_matches[m_index];
            const double & focal1 = camera_params[m1].focal;
            const double & focal2 = camera_params[m2].focal;
            
            MatrixXd A = MatrixXd::Zero(matches_info.num_inliers * DIMENSION_2D,
                                        HOMOGRAPHY_VARIABLES_COUNT);
            
            for(int j = 0; j < matches_info.num_inliers; ++j) {
                Point2d p1 = Point2d(images_features[m1].keypoints[matches_info.matches[j].queryIdx].pt) -
                             Point2d(translation_matrix[m1].at<double>(0, 2), translation_matrix[m1].at<double>(1, 2));
                Point2d p2 = Point2d(images_features[m2].keypoints[matches_info.matches[j].trainIdx].pt) -
                             Point2d(translation_matrix[m2].at<double>(0, 2), translation_matrix[m2].at<double>(1, 2));
                A(2*j  , 0) =  p1.x;
                A(2*j  , 1) =  p1.y;
                A(2*j  , 2) =         focal1;
                A(2*j  , 6) = -p2.x *   p1.x / focal2;
                A(2*j  , 7) = -p2.x *   p1.y / focal2;
                A(2*j  , 8) = -p2.x * focal1 / focal2;
                
                A(2*j+1, 3) =  p1.x;
                A(2*j+1, 4) =  p1.y;
                A(2*j+1, 5) =         focal1;
                A(2*j+1, 6) = -p2.y *   p1.x / focal2;
                A(2*j+1, 7) = -p2.y *   p1.y / focal2;
                A(2*j+1, 8) = -p2.y * focal1 / focal2;
            }
            JacobiSVD<MatrixXd, HouseholderQRPreconditioner> jacobi_svd(A, ComputeThinV);
            MatrixXd V = jacobi_svd.matrixV();
            Mat R(3, 3, CV_64FC1);
            for(int j = 0; j < V.rows(); ++j) {
                R.at<double>(j / 3, j % 3) = V(j, V.rows() - 1);
            }
            SVD svd(R, SVD::FULL_UV);
            relative_3D_rotations[m1][m2] = svd.u * svd.vt;
        }
        queue<int> que;
        vector<bool> labels(images_data.size(), false);
        const int & center_index = parameter.center_image_index;
        const vector<vector<bool> > & images_match_graph = parameter.getImagesMatchGraph();
        
        que.push(center_index);
        relative_3D_rotations[center_index][center_index] = Mat::eye(3, 3, CV_64FC1);
        
        while(que.empty() == false) {
            int now = que.front();
            que.pop();
            labels[now] = true;
            for(int i = 0; i < images_data.size(); ++i) {
                if(labels[i] == false) {
                    if(images_match_graph[now][i]) {
                        relative_3D_rotations[i][i] = relative_3D_rotations[now][i] * relative_3D_rotations[now][now];
                        que.push(i);
                    }
                    if(images_match_graph[i][now]) {
                        relative_3D_rotations[i][i] = relative_3D_rotations[i][now].inv() * relative_3D_rotations[now][now];
                        que.push(i);
                    }
                }
            }
        }
        /********************/
        for(int i = 0; i < camera_params.size(); ++i) {
            camera_params[i].aspect = 1;
            camera_params[i].ppx = translation_matrix[i].at<double>(0, 2);
            camera_params[i].ppy = translation_matrix[i].at<double>(1, 2);
            camera_params[i].t = Mat::zeros(3, 1, CV_64FC1);
            camera_params[i].R = relative_3D_rotations[i][i].inv();
            camera_params[i].R.convertTo(camera_params[i].R, CV_32FC1);
        }
        
        Ptr<detail::BundleAdjusterBase> adjuster = makePtr<detail::BundleAdjusterReproj>();
        adjuster->setTermCriteria(TermCriteria(TermCriteria::EPS, CRITERIA_MAX_COUNT, CRITERIA_EPSILON));
        
        Mat_<uchar> refine_mask = Mat::zeros(3, 3, CV_8U);
        refine_mask(0, 0) = 1; /* (0, 0)->focal, (0, 2)->ppx, (1, 2)->ppy, (1, 1)->aspect */
        adjuster->setConfThresh(1.f);
        adjuster->setRefinementMask(refine_mask);
        
        if (!(*adjuster)(images_features, pairwise_matches, camera_params)) {
            printError("F(getCameraParams) camera parameters adjuster failed");
        }
        
        Mat center_rotation_inv = camera_params[parameter.center_image_index].R.inv();
        for(int i = 0; i < camera_params.size(); ++i) {
            camera_params[i].R = center_rotation_inv * camera_params[i].R;
        }
        /* wave correction */
        if(WAVE_CORRECT != WAVE_X) {
            vector<Mat> rotations;
            rotations.reserve(camera_params.size());
            for(int i = 0; i < camera_params.size(); ++i) {
                rotations.emplace_back(camera_params[i].R);
            }
            waveCorrect(rotations, ((WAVE_CORRECT == WAVE_H) ? detail::WAVE_CORRECT_HORIZ : detail::WAVE_CORRECT_VERT));
            for(int i = 0; i < camera_params.size(); ++i) {
                camera_params[i].R = rotations[i];
            }
        }
        /*******************/
    }
    return camera_params;
}

const vector<vector<bool> > & MultiImages::getImagesFeaturesMaskByMatchingPoints() const {
    if(images_features_mask.empty()) {
        doFeatureMatching();
    }
    return images_features_mask;
}

const vector<vector<vector<bool> > > & MultiImages::getAPAPOverlapMask() const {
    if(apap_overlap_mask.empty()) {
        doFeatureMatching();
    }
    return apap_overlap_mask;
}
const vector<vector<vector<Mat> > > & MultiImages::getAPAPHomographies() const {
    if(apap_homographies.empty()) {
        doFeatureMatching();
    }
    return apap_homographies;
}

const vector<vector<vector<Point2> > > & MultiImages::getAPAPMatchingPoints() const {
    if(apap_matching_points.empty()) {
        doFeatureMatching();
    }
    return apap_matching_points;
}

const vector<vector<InterpolateVertex> > & MultiImages::getInterpolateVerticesOfMatchingPoints() const {
    if(mesh_interpolate_vertex_of_matching_pts.empty()) {
        mesh_interpolate_vertex_of_matching_pts.resize(images_data.size());
        const vector<detail::ImageFeatures> & images_features = getImagesFeaturesByMatchingPoints();
        for(int i = 0; i < mesh_interpolate_vertex_of_matching_pts.size(); ++i) {
            mesh_interpolate_vertex_of_matching_pts[i].reserve(images_features[i].keypoints.size());
            for(int j = 0; j < images_features[i].keypoints.size(); ++j) {
                mesh_interpolate_vertex_of_matching_pts[i].emplace_back(images_data[i].mesh_2d->getInterpolateVertex(images_features[i].keypoints[j].pt));
            }
        }
    }
    return mesh_interpolate_vertex_of_matching_pts;
}

const vector<int> & MultiImages::getImagesVerticesStartIndex() const {
    if(images_vertices_start_index.empty()) {
        images_vertices_start_index.reserve(images_data.size());
        int index = 0;
        for(int i = 0; i < images_data.size(); ++i) {
            images_vertices_start_index.emplace_back(index);
            index += images_data[i].mesh_2d->getVertices().size() * DIMENSION_2D;
        }
    }
    return images_vertices_start_index;
}

const vector<SimilarityElements> & MultiImages::getImagesSimilarityElements(const enum GLOBAL_ROTATION_METHODS & _global_rotation_method) const {
    const vector<vector<SimilarityElements> *> & images_similarity_elements = {
        &images_similarity_elements_2D, &images_similarity_elements_3D
    };
    vector<SimilarityElements> & result = *images_similarity_elements[_global_rotation_method];
    if(result.empty()) {
        result.reserve(images_data.size());
        const vector<detail::CameraParams> & camera_params = getCameraParams();
        for(int i = 0; i < images_data.size(); ++i) {
            result.emplace_back(fabs(camera_params[parameter.center_image_index].focal / camera_params[i].focal),
                                -getEulerZXYRadians<float>(camera_params[i].R)[2]);
        }
        double rotate_theta = parameter.center_image_rotation_angle;
        for(int i = 0; i < images_data.size(); ++i) {
            double a = (result[i].theta - rotate_theta) * 180 / M_PI;
            result[i].theta = normalizeAngle(a) * M_PI / 180;
        }
        
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        const vector<vector<pair<double, double> > > & images_relative_rotation_range = getImagesRelativeRotationRange();

        switch (_global_rotation_method) {
            case GLOBAL_ROTATION_2D_METHOD:
            {
                class RotationNode {
                public:
                    int index, parent;
                    RotationNode(const int _index, const int _parent) {
                        index = _index, parent = _parent;
                    }
                private:
                    
                };
                const double TOLERANT_THETA = TOLERANT_ANGLE * M_PI / 180;
                vector<pair<int, double> > theta_constraints;
                vector<bool> decided(images_data.size(), false);
                vector<RotationNode> priority_que;
                theta_constraints.emplace_back(parameter.center_image_index, result[parameter.center_image_index].theta);
                decided[parameter.center_image_index] = true;
                priority_que.emplace_back(parameter.center_image_index, -1);
                const vector<vector<bool> > & images_match_graph = parameter.getImagesMatchGraph();
                while(priority_que.empty() == false) {
                    RotationNode node = priority_que.front();
                    priority_que.erase(priority_que.begin());
                    if(!decided[node.index]) {
                        decided[node.index] = true;
                        result[node.index].theta = result[node.parent].theta + getImagesMinimumLineDistortionRotation(node.parent, node.index);
                    }
                    for(int i = 0; i < decided.size(); ++i) {
                        if(!decided[i]) {
                            const int e[EDGE_VERTEX_SIZE] = { node.index, i };
                            for(int j = 0; j < EDGE_VERTEX_SIZE; ++j) {
                                if(images_match_graph[e[j]][e[!j]]) {
                                    RotationNode new_node(i, node.index);
                                    if(isRotationInTheRange<double>(0, result[node.index].theta + images_relative_rotation_range[node.index][i].first  - TOLERANT_THETA,
                                                                       result[node.index].theta + images_relative_rotation_range[node.index][i].second + TOLERANT_THETA)) {
                                        priority_que.insert(priority_que.begin(), new_node);
                                        result[i].theta = 0;
                                        decided[i] = true;
                                        theta_constraints.emplace_back(i, 0);
                                    } else {
                                        priority_que.emplace_back(new_node);
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
                const int equations_count = (int)(images_match_graph_pair_list.size() + theta_constraints.size()) * DIMENSION_2D;
                SparseMatrix<double> A(equations_count, images_data.size() * DIMENSION_2D);
                VectorXd b = VectorXd::Zero(equations_count);
                vector<Triplet<double> > triplets;
                triplets.reserve(theta_constraints.size() * 2 + images_match_graph_pair_list.size() * 6);
                
                int equation = 0;
                for(int i = 0; i < theta_constraints.size(); ++i) {
                    triplets.emplace_back(equation    , DIMENSION_2D * theta_constraints[i].first    , STRONG_CONSTRAINT);
                    triplets.emplace_back(equation + 1, DIMENSION_2D * theta_constraints[i].first + 1, STRONG_CONSTRAINT);
                    b[equation    ] = STRONG_CONSTRAINT * cos(theta_constraints[i].second);
                    b[equation + 1] = STRONG_CONSTRAINT * sin(theta_constraints[i].second);
                    equation += DIMENSION_2D;
                } 
                for(int i = 0; i < images_match_graph_pair_list.size(); ++i) {
                    const pair<int, int> & match_pair = images_match_graph_pair_list[i];
                    const int & m1 = match_pair.first, & m2 = match_pair.second;
                    const FLOAT_TYPE & MLDR_theta = getImagesMinimumLineDistortionRotation(m1, m2);
                    triplets.emplace_back(equation    , DIMENSION_2D * m1    ,  cos(MLDR_theta));
                    triplets.emplace_back(equation    , DIMENSION_2D * m1 + 1, -sin(MLDR_theta));
                    triplets.emplace_back(equation    , DIMENSION_2D * m2    ,               -1);
                    triplets.emplace_back(equation + 1, DIMENSION_2D * m1    ,  sin(MLDR_theta));
                    triplets.emplace_back(equation + 1, DIMENSION_2D * m1 + 1,  cos(MLDR_theta));
                    triplets.emplace_back(equation + 1, DIMENSION_2D * m2 + 1,               -1);
                    equation += DIMENSION_2D;
                }
                assert(equation == equations_count);
                A.setFromTriplets(triplets.begin(), triplets.end());
                LeastSquaresConjugateGradient<SparseMatrix<double> > lscg(A);
                VectorXd x = lscg.solve(b);

               for(int i = 0; i < images_data.size(); ++i) {
                    result[i].theta = atan2(x[DIMENSION_2D * i + 1], x[DIMENSION_2D * i]);
                }
            }
                break;
            case GLOBAL_ROTATION_3D_METHOD:
            {
                const int equations_count = (int)images_match_graph_pair_list.size() * DIMENSION_2D + DIMENSION_2D;
                SparseMatrix<double> A(equations_count, images_data.size() * DIMENSION_2D);
                VectorXd b = VectorXd::Zero(equations_count);
                vector<Triplet<double> > triplets;
                triplets.reserve(images_match_graph_pair_list.size() * 6 + DIMENSION_2D);
                
                b[0] = STRONG_CONSTRAINT * cos(result[parameter.center_image_index].theta);
                b[1] = STRONG_CONSTRAINT * sin(result[parameter.center_image_index].theta);
                triplets.emplace_back(0, DIMENSION_2D * parameter.center_image_index    , STRONG_CONSTRAINT);
                triplets.emplace_back(1, DIMENSION_2D * parameter.center_image_index + 1, STRONG_CONSTRAINT);
                int equation = DIMENSION_2D;
                for(int i = 0; i < images_match_graph_pair_list.size(); ++i) {
                    const pair<int, int> & match_pair = images_match_graph_pair_list[i];
                    const int & m1 = match_pair.first, & m2 = match_pair.second;
                    const double guess_theta = result[m2].theta - result[m1].theta;
                    FLOAT_TYPE decision_theta, weight;
                    if(isRotationInTheRange(guess_theta,
                                            images_relative_rotation_range[m1][m2].first,
                                            images_relative_rotation_range[m1][m2].second)) {
                        decision_theta = guess_theta;
                        weight = LAMBDA_GAMMA;
                    } else {
                        decision_theta = getImagesMinimumLineDistortionRotation(m1, m2);
                        weight = 1;
                    }
                    triplets.emplace_back(equation    , DIMENSION_2D * m1    , weight *  cos(decision_theta));
                    triplets.emplace_back(equation    , DIMENSION_2D * m1 + 1, weight * -sin(decision_theta));
                    triplets.emplace_back(equation    , DIMENSION_2D * m2    ,                       -weight);
                    triplets.emplace_back(equation + 1, DIMENSION_2D * m1    , weight *  sin(decision_theta));
                    triplets.emplace_back(equation + 1, DIMENSION_2D * m1 + 1, weight *  cos(decision_theta));
                    triplets.emplace_back(equation + 1, DIMENSION_2D * m2 + 1,                       -weight);
                    
                    equation += DIMENSION_2D;
                }
                assert(equation == equations_count);
                A.setFromTriplets(triplets.begin(), triplets.end());
                LeastSquaresConjugateGradient<SparseMatrix<double> > lscg(A);
                VectorXd x = lscg.solve(b);
                
                for(int i = 0; i < images_data.size(); ++i) {
                    result[i].theta = atan2(x[DIMENSION_2D * i + 1], x[DIMENSION_2D * i]);
                }
            }
                break;
            default:
                printError("F(getImagesSimilarityElements) NISwGSP_ROTATION_METHOD");
                break;
        }
    }
    return result;
}

const vector<vector<pair<double, double> > > & MultiImages::getImagesRelativeRotationRange() const {
    if(images_relative_rotation_range.empty()) {
        images_relative_rotation_range.resize(images_data.size());
        for(int i = 0; i < images_relative_rotation_range.size(); ++i) {
            images_relative_rotation_range[i].resize(images_relative_rotation_range.size(), make_pair(0, 0));
        }
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        const vector<vector<vector<bool> > > & apap_overlap_mask = getAPAPOverlapMask();
        const vector<vector<vector<Point2> > > & apap_matching_points = getAPAPMatchingPoints();
        for(int i = 0; i < images_match_graph_pair_list.size(); ++i) {
            const pair<int, int> & match_pair = images_match_graph_pair_list[i];
            const int & m1 = match_pair.first, & m2 = match_pair.second;
            const vector<Edge> & m1_edges = images_data[m1].mesh_2d->getEdges();
            const vector<Edge> & m2_edges = images_data[m2].mesh_2d->getEdges();
            const vector<const vector<Edge> *> & edges = { &m1_edges, &m2_edges };
            const vector<pair<int, int> > pair_index = { make_pair(m1, m2), make_pair(m2, m1) };
            const vector<pair<const vector<Point2> *, const vector<Point2> *> > & vertices_pair = {
                make_pair(&images_data[m1].mesh_2d->getVertices(), &apap_matching_points[m1][m2]),
                make_pair(&images_data[m2].mesh_2d->getVertices(), &apap_matching_points[m2][m1])
            };
            vector<double> positive, negative;
            const vector<bool> sign_mapping = { false, true, true, false };
            for(int j = 0; j < edges.size(); ++j) {
                for(int k = 0; k < edges[j]->size(); ++k) {
                    const Edge & e = (*edges[j])[k];
                    if(apap_overlap_mask[pair_index[j].first][pair_index[j].second][e.indices[0]] &&
                       apap_overlap_mask[pair_index[j].first][pair_index[j].second][e.indices[1]]) {
                        const Point2d a = (*vertices_pair[j].first )[e.indices[0]] - (*vertices_pair[j].first )[e.indices[1]];
                        const Point2d b = (*vertices_pair[j].second)[e.indices[0]] - (*vertices_pair[j].second)[e.indices[1]];
                        const double theta = acos(a.dot(b) / (norm(a) * norm(b)));
                        const double direction = a.x * b.y - a.y * b.x;
                        int map = ((direction > 0) << 1) + j;
                        if(sign_mapping[map]) {
                            positive.emplace_back( theta);
                        } else {
                            negative.emplace_back(-theta);
                        }
                    }
                }
            }
            std::sort(positive.begin(), positive.end());
            std::sort(negative.begin(), negative.end());
            
            if(positive.empty() == false && negative.empty() == false) {
                if(positive.back() - negative.front() < M_PI) {
                    images_relative_rotation_range[m1][m2].first  = negative.front() + 2 * M_PI;
                    images_relative_rotation_range[m1][m2].second = positive.back()  + 2 * M_PI;
                    images_relative_rotation_range[m2][m1].first  = 2 * M_PI - positive.back();
                    images_relative_rotation_range[m2][m1].second = 2 * M_PI - negative.front();
                } else {
                    images_relative_rotation_range[m1][m2].first  = positive.front();
                    images_relative_rotation_range[m1][m2].second = negative.back()  + 2 * M_PI;
                    images_relative_rotation_range[m2][m1].first  =          - negative.back();
                    images_relative_rotation_range[m2][m1].second = 2 * M_PI - positive.front();

                }
            } else if(positive.empty() == false) {
                images_relative_rotation_range[m1][m2].first  =            positive.front();
                images_relative_rotation_range[m1][m2].second =            positive.back();
                images_relative_rotation_range[m2][m1].first  = 2 * M_PI - positive.back();
                images_relative_rotation_range[m2][m1].second = 2 * M_PI - positive.front();
            } else {
                images_relative_rotation_range[m1][m2].first  =  negative.front() + 2 * M_PI;
                images_relative_rotation_range[m1][m2].second =  negative.back()  + 2 * M_PI;
                images_relative_rotation_range[m2][m1].first  = -negative.back();
                images_relative_rotation_range[m2][m1].second = -negative.front();
            }
        }
    }
    return images_relative_rotation_range;
}

FLOAT_TYPE MultiImages::getImagesMinimumLineDistortionRotation(const int _from, const int _to) const {
    if(images_minimum_line_distortion_rotation.empty()) {
        images_minimum_line_distortion_rotation.resize(images_data.size());
        for(int i = 0; i < images_minimum_line_distortion_rotation.size(); ++i) {
            images_minimum_line_distortion_rotation[i].resize(images_data.size(), std::numeric_limits<float>::max());
        }
    }
    if(images_minimum_line_distortion_rotation[_from][_to] == std::numeric_limits<float>::max()) {
        const vector<LineData> & from_lines   = images_data[_from].getLines();
        const vector<LineData> &   to_lines   = images_data[_to  ].getLines();
        const vector<Point2>   & from_project = getImagesLinesProject(_from, _to);
        const vector<Point2>   &   to_project = getImagesLinesProject(_to, _from);
        
        const vector<const vector<LineData> *> & lines    = { &from_lines,   &to_lines   };
        const vector<const vector<Point2  > *> & projects = { &from_project, &to_project };
        const vector<int> & img_indices = { _to, _from };
        const vector<int> sign_mapping = { -1, 1, 1, -1 };
        
        vector<pair<double, double> > theta_weight_pairs;
        for(int i = 0; i < lines.size(); ++i) {
            const int & rows = images_data[img_indices[i]].img.rows;
            const int & cols = images_data[img_indices[i]].img.cols;
            const vector<pair<Point2, Point2> > & boundary_edgs = {
                make_pair(Point2(0,       0), Point2(cols,    0)),
                make_pair(Point2(cols,    0), Point2(cols, rows)),
                make_pair(Point2(cols, rows), Point2(   0, rows)),
                make_pair(Point2(   0, rows), Point2(   0,    0))
            };
            for(int j = 0; j < lines[i]->size(); ++j) {
                const Point2 & p1 = (*projects[i])[EDGE_VERTEX_SIZE * j    ];
                const Point2 & p2 = (*projects[i])[EDGE_VERTEX_SIZE * j + 1];
                const bool p1_in_img = (p1.x >= 0 && p1.x <= cols && p1.y >= 0 && p1.y <= rows);
                const bool p2_in_img = (p2.x >= 0 && p2.x <= cols && p2.y >= 0 && p2.y <= rows);
                
                const bool p_in_img[EDGE_VERTEX_SIZE] = { p1_in_img, p2_in_img };
                
                Point2 p[EDGE_VERTEX_SIZE] = { p1, p2 };
            
                if(!p1_in_img || !p2_in_img) {
                    vector<double> scales;
                    for(int k = 0; k < boundary_edgs.size(); ++k) {
                        double s1;
                        if(isEdgeIntersection(p1, p2, boundary_edgs[k].first, boundary_edgs[k].second, &s1)) {
                            scales.emplace_back(s1);
                        }
                    }
                    assert(scales.size() <= EDGE_VERTEX_SIZE);
                    if(scales.size() == EDGE_VERTEX_SIZE) {
                        assert(!p1_in_img && !p2_in_img);
                        if(scales.front() > scales.back()) {
                            iter_swap(scales.begin(), scales.begin() + 1);
                        }
                        for(int k = 0; k < scales.size(); ++k) {
                            p[k] = p1 + scales[k] * (p2 - p1);
                        }
                    } else if(!scales.empty()){
                        for(int k = 0; k < EDGE_VERTEX_SIZE; ++k) {
                            if(!p_in_img[k]) {
                                p[k] = p1 + scales.front() * (p2 - p1);
                            }
                        }
                    } else {
                        continue;
                    }
                }
                const Point2d a = (*lines[i])[j].data[1] - (*lines[i])[j].data[0];
                const Point2d b = p2 - p1;
                const double theta = acos(a.dot(b) / (norm(a) * norm(b)));
                const double direction = a.x * b.y - a.y * b.x;
                const int map = ((direction > 0) << 1) + i;
                const double b_length_2 = sqrt(b.x * b.x + b.y * b.y);
                theta_weight_pairs.emplace_back(theta * sign_mapping[map],
                                                (*lines[i])[j].length * (*lines[i])[j].width * b_length_2);
            }
        }
        Point2 dir(0, 0);
        for(int i = 0; i < theta_weight_pairs.size(); ++i) {
            const double & theta = theta_weight_pairs[i].first;
            dir += (theta_weight_pairs[i].second * Point2(cos(theta), sin(theta)));
        }
        images_minimum_line_distortion_rotation[_from][_to] = acos(dir.x / (norm(dir))) * (dir.y > 0 ? 1 : -1);
        images_minimum_line_distortion_rotation[_to][_from] = -images_minimum_line_distortion_rotation[_from][_to];
    }
    return images_minimum_line_distortion_rotation[_from][_to];
}

const vector<Point2> & MultiImages::getImagesLinesProject(const int _from, const int _to) const {
    if(images_lines_projects.empty()) {
        images_lines_projects.resize(images_data.size());
        for(int i = 0; i < images_lines_projects.size(); ++i) {
            images_lines_projects[i].resize(images_data.size());
        }
    }
    if(images_lines_projects[_from][_to].empty()) {
        const vector<vector<vector<Point2> > > & feature_matches = getFeatureMatches();
        const vector<LineData> & lines = images_data[_from].getLines();
        vector<Point2> points, project_points;
        points.reserve(lines.size() * EDGE_VERTEX_SIZE);
        for(int i = 0; i < lines.size(); ++i) {
            for(int j = 0; j < EDGE_VERTEX_SIZE; ++j) {
                points.emplace_back(lines[i].data[j]);
            }
        }
        vector<Mat> not_be_used;
        APAP_Stitching::apap_project(feature_matches[_from][_to], feature_matches[_to][_from], points, images_lines_projects[_from][_to], not_be_used);
    }
    return images_lines_projects[_from][_to];
}

const vector<Mat> & MultiImages::getImages() const {
    if(images.empty()) {
        images.reserve(images_data.size());
        for(int i = 0; i < images_data.size(); ++i) {
            images.emplace_back(images_data[i].img);
        }
    }
    return images;
}

class dijkstraNode {
public:
    int from, pos;
    double dis;
    dijkstraNode(const int & _from,
                 const int & _pos,
                 const double & _dis) : from(_from), pos(_pos), dis(_dis) {
    }
    bool operator < (const dijkstraNode & rhs) const {
        return dis > rhs.dis;
    }
};
const vector<vector<double> > & MultiImages::getImagesGridSpaceMatchingPointsWeight(const double _global_weight_gamma) const {
    if(_global_weight_gamma && images_polygon_space_matching_pts_weight.empty()) {
        images_polygon_space_matching_pts_weight.resize(images_data.size());
        const vector<vector<bool > > & images_features_mask = getImagesFeaturesMaskByMatchingPoints();
        const vector<vector<InterpolateVertex> > & mesh_interpolate_vertex_of_matching_pts = getInterpolateVerticesOfMatchingPoints();
        for(int i = 0; i < images_polygon_space_matching_pts_weight.size(); ++i) {
            const int polygons_count = (int)images_data[i].mesh_2d->getPolygonsIndices().size();
            vector<bool> polygons_has_matching_pts(polygons_count, false);
            for(int j = 0; j < images_features_mask[i].size(); ++j) {
                if(images_features_mask[i][j]) {
                    polygons_has_matching_pts[mesh_interpolate_vertex_of_matching_pts[i][j].polygon] = true;
                }
            }
            images_polygon_space_matching_pts_weight[i].reserve(polygons_count);
            priority_queue<dijkstraNode> que;
            
            for(int j = 0; j < polygons_has_matching_pts.size(); ++j) {
                if(polygons_has_matching_pts[j]) {
                    polygons_has_matching_pts[j] = false;
                    images_polygon_space_matching_pts_weight[i].emplace_back(0);
                    que.push(dijkstraNode(j, j, 0));
                } else {
                    images_polygon_space_matching_pts_weight[i].emplace_back(MAXFLOAT);
                }
            }
            const vector<Indices> & polygons_neighbors = images_data[i].mesh_2d->getPolygonsNeighbors();
            const vector<Point2> & polygons_center = images_data[i].mesh_2d->getPolygonsCenter();
            while(que.empty() == false) {
                const dijkstraNode now = que.top();
                const int index = now.pos;
                que.pop();
                if(polygons_has_matching_pts[index] == false) {
                    polygons_has_matching_pts[index] = true;
                    for(int j = 0; j < polygons_neighbors[index].indices.size(); ++j) {
                        const int n = polygons_neighbors[index].indices[j];
                        if(polygons_has_matching_pts[n] == false) {
                            const double dis = norm(polygons_center[n] - polygons_center[now.from]);
                            if(images_polygon_space_matching_pts_weight[i][n] > dis) {
                                images_polygon_space_matching_pts_weight[i][n] = dis;
                                que.push(dijkstraNode(now.from, n, dis));
                            }
                        }
                    }
                }
            }
            const double normalize_inv = 1. / norm(Point2i(images_data[i].img.cols, images_data[i].img.rows));
            for(int j = 0; j < images_polygon_space_matching_pts_weight[i].size(); ++j) {
                images_polygon_space_matching_pts_weight[i][j] = images_polygon_space_matching_pts_weight[i][j] * normalize_inv;
            }
        }
    }
    return images_polygon_space_matching_pts_weight;
}

void MultiImages::initialFeaturePairsSpace() const {
    feature_pairs.resize(images_data.size());
    for(int i = 0; i < images_data.size(); ++i) {
        feature_pairs[i].resize(images_data.size());
    }
}

void MultiImages::initialLineFeatureMatchesSpace() const {
    line_feature_pairs.resize(images_data.size());
    for(int i = 0; i < images_data.size(); ++i) {
        line_feature_pairs[i].resize(images_data.size());
    }
}

void MultiImages::initialOverlapRect() const {
		overlap_rect.resize(images_data.size());
		for (int i = 0; i < images_data.size(); i++)
			overlap_rect[i].resize(images_data.size());
}

const vector<vector<vector<pair<int, int> > > > & MultiImages::getFeaturePairs() const {
    if(feature_pairs.empty()) {
        initialFeaturePairsSpace();
        initialOverlapRect();
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        for(int i = 0; i < images_match_graph_pair_list.size(); ++i) {
            const pair<int, int> & match_pair = images_match_graph_pair_list[i];
            const vector<pair<int, int> > & initial_indices = getInitialFeaturePairs(match_pair);
            const vector<Point2> & m1_fpts = images_data[match_pair.first ].getFeaturePoints();
            const vector<Point2> & m2_fpts = images_data[match_pair.second].getFeaturePoints();
            vector<Point2> X, Y;
            X.reserve(initial_indices.size());
            Y.reserve(initial_indices.size());
            for(int j = 0; j < initial_indices.size(); ++j) {
                const pair<int, int> it = initial_indices[j];
                X.emplace_back(m1_fpts[it.first ]);
                Y.emplace_back(m2_fpts[it.second]);
            }
            vector<pair<int, int> > & result = feature_pairs[match_pair.first][match_pair.second];
            Mat transform_matrix;
            result = getFeaturePairsBySequentialRANSAC(match_pair, X, Y, initial_indices, transform_matrix);
            double tx = transform_matrix.at<double>(0, 2);
			double ty = transform_matrix.at<double>(1, 2);
			Point2d corner1 = Point2d(tx, ty);
			Point2d corner2 = Point2d(0.0, 0.0);
			double x_tl = std::max(corner1.x, corner2.x);
			double y_tl = std::max(corner1.y, corner2.y);
			double x_br = std::min(corner1.x + images_data[match_pair.first].img.cols, corner2.x + images_data[match_pair.second].img.cols);
			double y_br = std::min(corner1.y + images_data[match_pair.first].img.rows, corner2.y + images_data[match_pair.second].img.rows);
            //计算的是第二幅图像到第一幅图像的H矩阵，据此可以计算出重叠区域的坐标
			//重叠区域在第一幅图像中的位置
			overlap_rect[match_pair.first][match_pair.second] = Rect2d(x_tl - tx, y_tl - ty, x_br - x_tl, y_br - y_tl);
			//重叠区域在第二幅图像中的位置
			overlap_rect[match_pair.second][match_pair.first] = Rect2d(x_tl, y_tl, x_br - x_tl, y_br - y_tl);
            assert(result.empty() == false);
#ifndef NDEBUG
            writeImageOfFeaturePairs("sRANSAC", match_pair, result);
            // test line for line feature RANSAC
            //const vector<vector<vector<DMatch> > > & LineFeatureRANSAC = getLineFeaturePairs();
#endif
        }
    }
    return feature_pairs;
}

//Added by @Bingkun, 30/3/2020
const vector<vector<vector<DMatch> > > & MultiImages::getLineFeaturePairs() const{
    if(line_feature_pairs.empty()){
        initialLineFeatureMatchesSpace();
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        for(int i = 0; i < images_match_graph_pair_list.size(); ++i){
            const pair<int, int> & match_pair = images_match_graph_pair_list[i];
            const vector<DMatch> & initialLineFeaturePairs = getInitialLineFeatureMatches(match_pair);
            const vector<KeyLine> & keylines1 = images_data[match_pair.first ].getKeyLines();
            const vector<KeyLine> & keylines2 = images_data[match_pair.second ].getKeyLines();
            vector<Point2f> X, Y;
            for(int j = 0; j < initialLineFeaturePairs.size(); ++j) {
                const DMatch it = initialLineFeaturePairs[j];
                X.emplace_back(keylines1[it.queryIdx].getStartPoint());
                Y.emplace_back(keylines2[it.trainIdx].getStartPoint());
            }
            vector<DMatch> & result = line_feature_pairs[match_pair.first][match_pair.second];
            result = getLineFeatureMatchesBySequentialRANSAC(match_pair, X, Y, initialLineFeaturePairs);
            assert(result.empty() == false);
            std::cout << "RANSAC_First Stage Line Binary Descriptor Matcher is : " << result.size() << std::endl;
            X.clear();
            Y.clear();
            for(int j = 0; j < result.size(); ++j) {
                const DMatch it = result[j];
                X.emplace_back(keylines1[it.queryIdx].getEndPoint());
                Y.emplace_back(keylines2[it.trainIdx].getEndPoint());
            }
            result = getLineFeatureMatchesBySequentialRANSAC(match_pair, X, Y, result);
            std::cout << "Line RANSAC Done !" << std::endl;

            Mat lineFeaturePairsOutImg;
	        vector<char> mask(matches.size(), 1);
	        drawLineMatches(images_data[match_pair.first ].img, lbd_octave1, images_data[match_pair.second ].img, lbd_octave2, result, lineFeaturePairsOutImg,
		            Scalar::all(-1), Scalar::all(-1), mask, DrawLinesMatchesFlags::DEFAULT);
	        std::cout << "RANSAC_Second Stage Line Binary Descriptor Matcher is : " << result.size() << std::endl;

            string _name = "RANSAC";
            imwrite(parameter.debug_dir +
            "LSD_line_feature_pairs-" + _name + "-" +
            images_data[match_pair.first ].file_name  + "-" +
            images_data[match_pair.second].file_name  + "-" +
            to_string(result.size()) +
            images_data[match_pair.first ].file_extension, lineFeaturePairsOutImg);
        }
    }
    return line_feature_pairs;
}

const vector<vector<vector<Point2> > > & MultiImages::getFeatureMatches() const {
    if(feature_matches.empty()) {
        const vector<vector<vector<pair<int, int> > > > & feature_pairs = getFeaturePairs();
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        feature_matches.resize(images_data.size());
        for(int i = 0; i < images_data.size(); ++i) {
            feature_matches[i].resize(images_data.size());
        }
        for(int i = 0; i < images_match_graph_pair_list.size(); ++i) {
            const pair<int, int> & match_pair = images_match_graph_pair_list[i];
            const int & m1 = match_pair.first, & m2 = match_pair.second;
            feature_matches[m1][m2].reserve(feature_pairs[m1][m2].size());
            feature_matches[m2][m1].reserve(feature_pairs[m1][m2].size());
            const vector<Point2> & m1_fpts = images_data[m1].getFeaturePoints();
            const vector<Point2> & m2_fpts = images_data[m2].getFeaturePoints();
            for(int j = 0; j < feature_pairs[m1][m2].size(); ++j) {
                feature_matches[m1][m2].emplace_back(m1_fpts[feature_pairs[m1][m2][j].first ]);
                feature_matches[m2][m1].emplace_back(m2_fpts[feature_pairs[m1][m2][j].second]);
            }
        }
    }
    return feature_matches;
}

//Added by @Bingkun, 14/5/2020
/*
 *手动选取直线
 */
vector<Point2f> temp_pt;
pair<Point2f, Point2f> lines_pt;
void onmouse1(int event, int x, int y, int flags, void *utsc){
    Mat &image = *(Mat*)utsc;

    if(event == CV_EVENT_LBUTTONDOWN){
        if(temp_pt.size() == 2){
            line(image,temp_pt[0],temp_pt[1],Scalar(0,0,255),1,LINE_AA);
            lines_pt = make_pair(temp_pt[0], temp_pt[1]);
            temp_pt.clear();
            return;
        }
        else{
            circle(image,Point2f(x,y),1,Scalar(0,0,255),2);
            temp_pt.emplace_back(Point2f(x, y));
            imshow("Lines1", image);
        }
    }
    else if(event == CV_EVENT_MBUTTONDOWN)
        temp_pt.clear();
}

void onmouse2(int event, int x, int y, int flags, void *utsc){
    Mat &image = *(Mat*)utsc;
    if(event == CV_EVENT_LBUTTONDOWN){
        if(temp_pt.size() == 2){
            line(image,temp_pt[0],temp_pt[1],Scalar(0,0,255),1,LINE_AA);
            lines_pt = make_pair(temp_pt[0], temp_pt[1]);
            temp_pt.clear();
            return;
        }
        else{
            circle(image,Point2f(x,y),1,Scalar(0,0,255),2);
            temp_pt.emplace_back(Point2f(x, y));
            imshow("Lines2", image);
        } 
    }
    else if(event == CV_EVENT_MBUTTONDOWN)
        temp_pt.clear();
}

const vector<vector<vector<KeyLine> > > & MultiImages::getLineFeatureMatches() const {
    if(line_feature_matches.empty()){
        const vector<vector<vector<DMatch> > > & line_feature_pairs = getLineFeaturePairs();
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        line_feature_matches.resize(images_data.size());
        for(int i = 0; i < images_data.size(); ++i) {
           line_feature_matches[i].resize(images_data.size());
        }
        for(int i = 0; i< images_match_graph_pair_list.size(); ++i){
            const pair<int, int> & match_pair = images_match_graph_pair_list[i];
            const int & m1 = match_pair.first, & m2 = match_pair.second;
            line_feature_matches[m1][m2].reserve(line_feature_pairs[m1][m2].size());
            line_feature_matches[m2][m1].reserve(line_feature_pairs[m1][m2].size());
            const vector<KeyLine> & keylines1 = images_data[m1].getKeyLines();
            const vector<KeyLine> & keylines2 = images_data[m2].getKeyLines();
            for(int j = 0; j < line_feature_pairs[m1][m2].size(); ++j){
                line_feature_matches[m1][m2].emplace_back(keylines1[line_feature_pairs[m1][m2][j].queryIdx]);
                line_feature_matches[m2][m1].emplace_back(keylines2[line_feature_pairs[m1][m2][j].trainIdx]);
            }

            Mat lines_image1 = images_data[m1].img.clone();
            Mat lines_image2 = images_data[m2].img.clone();
            cout << "Select Line Match manually? Press 'L'." << endl;
            string flag;
            getline(cin, flag);
            if(flag == "L" || flag == "l"){
                namedWindow("Lines1", 0);
                imshow("Lines1", lines_image1);
                cout << "Please select a line in the first pic, press any key after you done." << endl;
                setMouseCallback("Lines1", onmouse1,(void *) &lines_image1);
                waitKey();
                KeyLine manual_line1;
                manual_line1.startPointX = lines_pt.first.x;
                manual_line1.startPointY = lines_pt.first.y;
                manual_line1.endPointX = lines_pt.second.x;
                manual_line1.endPointY = lines_pt.second.y;
                line_feature_matches[m1][m2].emplace_back(manual_line1);
                images_data[m1].keylines.emplace_back(manual_line1);
                temp_pt.clear();
                imshow("Lines1", lines_image1);
                namedWindow("Lines2", 0);
                imshow("Lines2", lines_image2);
                cout << "Please select a line in the second pic, press any key after you done." << endl;
                setMouseCallback("Lines2", onmouse2,(void *) &lines_image2);
                waitKey();
                KeyLine manual_line2;
                manual_line2.startPointX = lines_pt.first.x;
                manual_line2.startPointY = lines_pt.first.y;
                manual_line2.endPointX = lines_pt.second.x;
                manual_line2.endPointY = lines_pt.second.y;
                line_feature_matches[m2][m1].emplace_back(manual_line2);
                images_data[m2].keylines.emplace_back(manual_line2);
                temp_pt.clear();
                destroyAllWindows();
            }
        }
    }
    return line_feature_matches;
}

vector<pair<int, int> > MultiImages::getFeaturePairsBySequentialRANSAC(const pair<int, int> & _match_pair,
                                                                       const vector<Point2> & _X,
                                                                       const vector<Point2> & _Y,
                                                                       const vector<pair<int, int> > & _initial_indices,
                                                                       Mat & _global_homography) const {
    const int HOMOGRAPHY_MODEL_MIN_POINTS = 4;
    const int GLOBAL_MAX_ITERATION = log(1 - OPENCV_DEFAULT_CONFIDENCE) / log(1 - pow(GLOBAL_TRUE_PROBABILITY, HOMOGRAPHY_MODEL_MIN_POINTS));
    
    vector<char> final_mask(_initial_indices.size(), 0);
    _global_homography = findHomography(_X, _Y, CV_RANSAC, parameter.global_homography_max_inliers_dist, final_mask, GLOBAL_MAX_ITERATION);
    
    vector<Point2> tmp_X = _X, tmp_Y = _Y;
    
    vector<int> mask_indices(_initial_indices.size(), 0);
    for(int i = 0; i < mask_indices.size(); ++i) {
        mask_indices[i] = i;
    }
    
    while(tmp_X.size() >= HOMOGRAPHY_MODEL_MIN_POINTS &&
          parameter.local_homogrpahy_max_inliers_dist < parameter.global_homography_max_inliers_dist) {
        const int LOCAL_MAX_ITERATION = log(1 - OPENCV_DEFAULT_CONFIDENCE) / log(1 - pow(LOCAL_TRUE_PROBABILITY, HOMOGRAPHY_MODEL_MIN_POINTS));
        vector<Point2> next_X, next_Y;
        vector<char> mask(tmp_X.size(), 0);
        findHomography(tmp_X, tmp_Y, CV_RANSAC, parameter.local_homogrpahy_max_inliers_dist, mask, LOCAL_MAX_ITERATION);

        int inliers_count = 0;
        for(int i = 0; i < mask.size(); ++i) {
            if(mask[i]) { ++inliers_count; }
        }
        if(inliers_count < parameter.local_homography_min_features_count) {
            break;
        }
        for(int i = 0, shift = -1; i < mask.size(); ++i) {
            if(mask[i]) {
                final_mask[mask_indices[i]] = 1;
            } else {
                next_X.emplace_back(tmp_X[i]);
                next_Y.emplace_back(tmp_Y[i]);
                mask_indices[++shift] = mask_indices[i];
            }
        }
#ifndef NDEBUG
        cout << "Local true Probabiltiy = " << next_X.size() / (float)tmp_X.size() << endl;
#endif
        tmp_X = next_X;
        tmp_Y = next_Y;
    }
    vector<pair<int, int> > result;
    for(int i = 0; i < final_mask.size(); ++i) {
        if(final_mask[i]) {
            result.emplace_back(_initial_indices[i]);
        }
    }
#ifndef NDEBUG
    cout << "Global true Probabiltiy = " << result.size() / (float)_initial_indices.size() << endl;
#endif
    return result;
}

bool compareFeaturePair(const FeatureDistance & fd_1, const FeatureDistance & fd_2) {
    return
    (fd_1.feature_index[0] == fd_2.feature_index[0]) ?
    (fd_1.feature_index[1]  < fd_2.feature_index[1]) :
    (fd_1.feature_index[0]  < fd_2.feature_index[0]) ;
}

//Added by @Bingkun, 31/3/2020
vector<DMatch> MultiImages::getLineFeatureMatchesBySequentialRANSAC(const pair<int, int> & _match_pair,
                                                            const vector<Point2f> & _X,
                                                            const vector<Point2f> & _Y,
                                                            const vector<DMatch> & _initialLineFeaturePairs) const{
    const int HOMOGRAPHY_MODEL_MIN_POINTS = 4;
    const int GLOBAL_MAX_ITERATION = log(1 - OPENCV_DEFAULT_CONFIDENCE) / log(1 - pow(GLOBAL_TRUE_PROBABILITY, HOMOGRAPHY_MODEL_MIN_POINTS));
    
    vector<char> final_mask(_initialLineFeaturePairs.size(), 0);
    findHomography(_X, _Y, CV_RANSAC, parameter.global_homography_max_inliers_dist, final_mask, GLOBAL_MAX_ITERATION);
    
    vector<Point2f> tmp_X = _X, tmp_Y = _Y;
    
    vector<int> mask_indices(_initialLineFeaturePairs.size(), 0);
    for(int i = 0; i < mask_indices.size(); ++i) {
        mask_indices[i] = i;
    }
    
    while(tmp_X.size() >= HOMOGRAPHY_MODEL_MIN_POINTS &&
          parameter.local_homogrpahy_max_inliers_dist < parameter.global_homography_max_inliers_dist) {
        const int LOCAL_MAX_ITERATION = log(1 - OPENCV_DEFAULT_CONFIDENCE) / log(1 - pow(LOCAL_TRUE_PROBABILITY, HOMOGRAPHY_MODEL_MIN_POINTS));
        vector<Point2f> next_X, next_Y;
        vector<char> mask(tmp_X.size(), 0);
        findHomography(tmp_X, tmp_Y, CV_RANSAC, parameter.local_homogrpahy_max_inliers_dist, mask, LOCAL_MAX_ITERATION);

        int inliers_count = 0;
        for(int i = 0; i < mask.size(); ++i) {
            if(mask[i]) { ++inliers_count; }
        }
        /*
        if(inliers_count < parameter.local_homography_min_features_count) {
            break;
        }
        */

        if(inliers_count < 15) {
            break;
        }
        for(int i = 0, shift = -1; i < mask.size(); ++i) {
            if(mask[i]) {
                final_mask[mask_indices[i]] = 1;
            } else {
                next_X.emplace_back(tmp_X[i]);
                next_Y.emplace_back(tmp_Y[i]);
                mask_indices[++shift] = mask_indices[i];
            }
        }
#ifndef NDEBUG
        cout << "Local true Probabiltiy = " << next_X.size() / (float)tmp_X.size() << endl;
#endif
        tmp_X = next_X;
        tmp_Y = next_Y;
    }
    vector<DMatch> result;
    for(int i = 0; i < final_mask.size(); ++i) {
        if(final_mask[i]) {
            result.emplace_back(_initialLineFeaturePairs[i]);
        }
    }
#ifndef NDEBUG
    cout << "Global true Probabiltiy = " << result.size() / (float)_initialLineFeaturePairs.size() << endl;
#endif
    return result;
}

vector<pair<int, int> > MultiImages::getInitialFeaturePairs(const pair<int, int> & _match_pair) const {
    const int nearest_size = 2, pair_count = 1;
    const bool ratio_test = true, intersect = true;
    
    assert(nearest_size > 0);
    
    const int feature_size_1 = (int)images_data[_match_pair.first ].getFeaturePoints().size();
    const int feature_size_2 = (int)images_data[_match_pair.second].getFeaturePoints().size();
    const int PAIR_COUNT = 2;
    const int feature_size[PAIR_COUNT] = { feature_size_1, feature_size_2 };
    const int pair_match[PAIR_COUNT] = { _match_pair.first , _match_pair.second };
    vector<FeatureDistance> feature_pairs[PAIR_COUNT];
    
    for(int p = 0; p < pair_count; ++p) {
        const int another_feature_size = feature_size[1 - p];
        const int nearest_k = min(nearest_size, another_feature_size);
        const vector<FeatureDescriptor> & feature_descriptors_1 = images_data[pair_match[ p]].getFeatureDescriptors();
        const vector<FeatureDescriptor> & feature_descriptors_2 = images_data[pair_match[!p]].getFeatureDescriptors();
        for(int f1 = 0; f1 < feature_size[p]; ++f1) {
            set<FeatureDistance> feature_distance_set;
            feature_distance_set.insert(FeatureDistance(MAXFLOAT, p, -1, -1));
            for(int f2 = 0; f2 < feature_size[!p]; ++f2) {
                const double dist = FeatureDescriptor::getDistance(feature_descriptors_1[f1], feature_descriptors_2[f2], feature_distance_set.begin()->distance);
                if(dist < feature_distance_set.begin()->distance) {
                    if(feature_distance_set.size() == nearest_k) {
                        feature_distance_set.erase(feature_distance_set.begin());
                    }
                    feature_distance_set.insert(FeatureDistance(dist, p, f1, f2));
                }
            }
            set<FeatureDistance>::const_iterator it = feature_distance_set.begin();
            if(ratio_test) {
                const set<FeatureDistance>::const_iterator it2 = std::next(it, 1);
                if(nearest_k == nearest_size &&
                   it2->distance * FEATURE_RATIO_TEST_THRESHOLD > it->distance) {
                    continue;
                }
                it = it2;
            }
            feature_pairs[p].insert(feature_pairs[p].end(), it, feature_distance_set.end());
        }
    }
    vector<FeatureDistance> feature_pairs_result;
    if(pair_count == PAIR_COUNT) {
        std::sort(feature_pairs[0].begin(), feature_pairs[0].end(), compareFeaturePair);
        std::sort(feature_pairs[1].begin(), feature_pairs[1].end(), compareFeaturePair);
        if(intersect) {
            set_intersection(feature_pairs[0].begin(), feature_pairs[0].end(),
                             feature_pairs[1].begin(), feature_pairs[1].end(),
                             std::inserter(feature_pairs_result, feature_pairs_result.begin()),
                             compareFeaturePair);
        } else {
            set_union(feature_pairs[0].begin(), feature_pairs[0].end(),
                      feature_pairs[1].begin(), feature_pairs[1].end(),
                      std::inserter(feature_pairs_result, feature_pairs_result.begin()),
                      compareFeaturePair);
        }
    } else {
        feature_pairs_result = std::move(feature_pairs[0]);
    }
    vector<double> distances;
    distances.reserve(feature_pairs_result.size());
    for(int i = 0; i < feature_pairs_result.size(); ++i) {
        distances.emplace_back(feature_pairs_result[i].distance);
    }
    double mean, std;
    Statistics::getMeanAndSTD(distances, mean, std);
    
    const double OUTLIER_THRESHOLD = (INLIER_TOLERANT_STD_DISTANCE * std) + mean;
    vector<pair<int, int> > initial_indices;
    initial_indices.reserve(feature_pairs_result.size());
    for(int i = 0; i < feature_pairs_result.size(); ++i) {
        if(feature_pairs_result[i].distance < OUTLIER_THRESHOLD) {
            initial_indices.emplace_back(feature_pairs_result[i].feature_index[0],
                                         feature_pairs_result[i].feature_index[1]);
        }
    }

#ifndef NDEBUG
    writeImageOfFeaturePairs("init", _match_pair, initial_indices);
    // test line for line feature extraction
    //vector<DMatch> line_feature_match = getInitialLineFeatureMatches(_match_pair);
#endif
    return initial_indices;
}

vector<DMatch> MultiImages::getInitialLineFeatureMatches(const pair<int, int> & _match_pair) const{
    const vector<KeyLine> & keylines1 = images_data[_match_pair.first ].getKeyLines();
    const vector<KeyLine> & keylines2 = images_data[_match_pair.second ].getKeyLines();
    const Mat & descr1 = images_data[_match_pair.first ].getLineFeatureDescriptor();
    const Mat & descr2 = images_data[_match_pair.second ].getLineFeatureDescriptor();
    
    Mat left_lbd, right_lbd;
    matches.clear();
    lbd_octave1.clear();
    lbd_octave2.clear();
    for (int i = 0; i < (int)keylines1.size(); i++)
	{
		if (keylines1[i].octave == 0)
		{
			lbd_octave1.push_back(keylines1[i]);
			left_lbd.push_back(descr1.row(i));
		}
	}

	for (int j = 0; j < (int)keylines2.size(); j++)
	{
		if (keylines2[j].octave == 0)
		{
			lbd_octave2.push_back(keylines2[j]);
			right_lbd.push_back(descr2.row(j));
		}
	}

    Ptr<BinaryDescriptorMatcher> bdm = BinaryDescriptorMatcher::createBinaryDescriptorMatcher();
    
	bdm->match(left_lbd, right_lbd, matches);
    vector<DMatch> good_matches;
	for (int i = 0; i < (int)matches.size(); i++)
	{
		if (matches[i].distance < MATCHES_DIST_THRESHOLD)
			good_matches.push_back(matches[i]);
	}
    
    Mat lineFeaturePairsOutImg;
	vector<char> mask(matches.size(), 1);
	drawLineMatches(images_data[_match_pair.first ].img, lbd_octave1, images_data[_match_pair.second ].img, lbd_octave2, good_matches, lineFeaturePairsOutImg,
		            Scalar::all(-1), Scalar::all(-1), mask, DrawLinesMatchesFlags::DEFAULT);
	cout << "Line Binary Descriptor Matcher is : " << good_matches.size() << endl;

    string _name = "init";
    imwrite(parameter.debug_dir +
            "LSD_line_feature_pairs-" + _name + "-" +
            images_data[_match_pair.first ].file_name  + "-" +
            images_data[_match_pair.second].file_name  + "-" +
            to_string(good_matches.size()) +
            images_data[_match_pair.first ].file_extension, lineFeaturePairsOutImg);

    return good_matches;
}

Mat MultiImages::textureMapping(const vector<vector<Point2> > & _vertices,
                                const Size2 & _target_size,
                                const BLENDING_METHODS & _blend_method) const {
    vector<Mat> warp_images;
    return textureMapping(_vertices, _target_size, _blend_method, warp_images);
}

Mat MultiImages::textureMapping(const vector<vector<Point2> > & _vertices,
                                const Size2 & _target_size,
                                const BLENDING_METHODS & _blend_method,
                                vector<Mat> & _warp_images) const {
    
    vector<Mat> weight_mask, new_weight_mask, warp_images_mask(images_data.size());
    vector<Point2> origins;
    vector<Rect_<FLOAT_TYPE> > rects = getVerticesRects<FLOAT_TYPE>(_vertices);
    
    switch (_blend_method) {
        case BLEND_AVERAGE:
            break;
        case BLEND_LINEAR:
            weight_mask = getMatsLinearBlendWeight(getImages());
            break;
        case BLEND_SEAM:
            break;
        default:
            printError("F(textureMapping) BLENDING METHOD");;
    }
#ifndef NDEBUG
    for(int i = 0; i < rects.size(); ++i) {
        cout << images_data[i].file_name << " rect = " << rects[i] << endl;
    }
#endif
    _warp_images.reserve(_vertices.size());
    origins.reserve(_vertices.size());
    new_weight_mask.reserve(_vertices.size());
    
    const int NO_GRID = -1, TRIANGLE_COUNT = 3, PRECISION = 0;
    const int SCALE = pow(2, PRECISION);
    
    for(int i = 0; i < images_data.size(); ++i) {
        const vector<Point2> & src_vertices = images_data[i].mesh_2d->getVertices();
        const vector<Indices> & polygons_indices = images_data[i].mesh_2d->getPolygonsIndices();
        const Point2 origin(rects[i].x, rects[i].y);
        const Point2 shift(0.5, 0.5);
        vector<Mat> affine_transforms;
        affine_transforms.reserve(polygons_indices.size() * (images_data[i].mesh_2d->getTriangulationIndices().size()));
        Mat polygon_index_mask(rects[i].height + shift.y, rects[i].width + shift.x, CV_32SC1, Scalar::all(NO_GRID));
        int label = 0;
        for(int j = 0; j < polygons_indices.size(); ++j) {
            for(int k = 0; k < images_data[i].mesh_2d->getTriangulationIndices().size(); ++k) {
                const Indices & index = images_data[i].mesh_2d->getTriangulationIndices()[k];
                const Point2i contour[] = {
                    (_vertices[i][polygons_indices[j].indices[index.indices[0]]] - origin) * SCALE,
                    (_vertices[i][polygons_indices[j].indices[index.indices[1]]] - origin) * SCALE,
                    (_vertices[i][polygons_indices[j].indices[index.indices[2]]] - origin) * SCALE,
                };
                fillConvexPoly(polygon_index_mask, contour, TRIANGLE_COUNT, label, LINE_AA, PRECISION);
                Point2f src[] = {
                    _vertices[i][polygons_indices[j].indices[index.indices[0]]] - origin,
                    _vertices[i][polygons_indices[j].indices[index.indices[1]]] - origin,
                    _vertices[i][polygons_indices[j].indices[index.indices[2]]] - origin
                };
                Point2f dst[] = {
                    src_vertices[polygons_indices[j].indices[index.indices[0]]],
                    src_vertices[polygons_indices[j].indices[index.indices[1]]],
                    src_vertices[polygons_indices[j].indices[index.indices[2]]]
                };
                affine_transforms.emplace_back(getAffineTransform(src, dst));
                ++label;
            }
        }
        Mat image = Mat::zeros(rects[i].height + shift.y, rects[i].width + shift.x, CV_8UC4);
        Mat mask = Mat::zeros(image.size(), CV_8UC1);//缝合线可以穿过的区域
        Mat w_mask = (_blend_method == BLEND_LINEAR) ? Mat::zeros(image.size(), CV_32FC1) : Mat();
        for(int y = 0; y < image.rows; ++y) {
            for(int x = 0; x < image.cols; ++x) {
                int polygon_index = polygon_index_mask.at<int>(y, x);
                if(polygon_index != NO_GRID) {
                    Point2 p_f = applyTransform2x3<FLOAT_TYPE>(x, y,
                                                               affine_transforms[polygon_index]);
                    if(p_f.x >= 0 && p_f.y >= 0 &&
                       p_f.x <= images_data[i].img.cols &&
                       p_f.y <= images_data[i].img.rows) {
                        Vec<uchar, 1> alpha = getSubpix<uchar, 1>(images_data[i].alpha_mask, p_f);
                        Vec3b c = getSubpix<uchar, 3>(images_data[i].img, p_f);
                        image.at<Vec4b>(y, x) = Vec4b(c[0], c[1], c[2], alpha[0]);
                        mask.at<uchar>(y, x) = 255;
                        if(_blend_method == BLEND_LINEAR) {
                            w_mask.at<float>(y, x) = getSubpix<float>(weight_mask[i], p_f);
                        }
                    }
                }
            }
        }
        /*
        polygon_index_mask.convertTo(polygon_index_mask, CV_8U);
		cvNamedWindow("polygon_index_mask", WINDOW_NORMAL);
		imshow("polygon_index_mask", polygon_index_mask);
		waitKey(1000);

        cout << "origin " << to_string(i + 1) << " : ( " << rects[i].x << ", " << rects[i].y << " )" << endl << endl;
		cout << "Img cols : " << images_data[i].img.cols << endl<<endl;
		namedWindow("_warp_images", WINDOW_NORMAL);
		imshow("_warp_images", image);
        waitKey(1000);
        */
        _warp_images.emplace_back(image);
        origins.emplace_back(rects[i].x, rects[i].y);
        warp_images_mask[i] = mask;
        if(_blend_method == BLEND_LINEAR) {
            new_weight_mask.emplace_back(w_mask);
        }
    }

    //cout << "TEST!!!!!!!!!!!"<<endl;
    std::vector<cv::UMat> U_warp_images(images_data.size()), U_new_weight_mask(images_data.size());
    std::vector<cv::Point> U_origins(images_data.size());
    for(int k = 0; k < U_warp_images.size(); k++){
        /*
        Mat img = Mat::zeros(_warp_images[k].size(), CV_32FC3);
        cvtColor(_warp_images[k], img, CV_RGBA2RGB);
        img.convertTo(img, CV_32FC3);
        img.copyTo(U_warp_images[k]);
        */
       _warp_images[k].copyTo(U_warp_images[k]);
    }
    for (int k = 0; k< U_new_weight_mask.size(); k++)
    {
        Mat warp_mask = Mat::zeros(warp_images_mask[k].size(), CV_8UC1);
        warp_images_mask[k].convertTo(warp_mask, CV_8UC1);
        warp_mask.copyTo(U_new_weight_mask[k]);
    }
    for (int k = 0; k< U_origins.size(); k++)
    {
        U_origins[k].x = floor(origins[k].x + 0.5);
        U_origins[k].y = floor(origins[k].y + 0.5);
    }

    Ptr<ExposureCompensator> compensator = 
                     ExposureCompensator::createDefault(ExposureCompensator::GAIN);//曝光补偿
    compensator->feed(U_origins, U_warp_images, U_new_weight_mask);
    for(int i=0; i<U_warp_images.size(); ++i)    //应用曝光补偿器，对图像进行曝光补偿
    {
        imwrite(parameter.debug_dir + "no_exposure_compensated" + to_string(i) + ".png", U_warp_images[i]);
        compensator->apply(i, U_origins[i], U_warp_images[i], U_new_weight_mask[i]);
        imwrite(parameter.debug_dir + "exposure_compensated" + to_string(i) + ".png", U_warp_images[i]);
        U_warp_images[i].copyTo(_warp_images[i]);
    }

    if(_blend_method == BLEND_SEAM){
        cout << "Seam-cutting Begins" <<endl;
        vector<UMat> U_warp_images_C3((int)(U_warp_images.size()));
        cout << "TEST!!!!!!" << endl;
        for(int k = 0; k < U_warp_images_C3.size(); k++){
            Mat img = Mat::zeros(U_warp_images[k].size(), CV_32FC3);
            cvtColor(U_warp_images[k], img, CV_RGBA2RGB);
            img.convertTo(img, CV_32FC3);
            img.copyTo(U_warp_images_C3[k]);
        }
        Ptr<SeamFinder> seam_finder = new GraphCutSeamFinder(GraphCutSeamFinderBase::COST_COLOR);
        seam_finder->find(U_warp_images_C3, U_origins, U_new_weight_mask);
        cout << "TEST2!!!!!!" << endl;
        int image_size = (int)(_warp_images[0].cols * _warp_images[0].rows);
        for(int k = 0; k < U_new_weight_mask.size(); k++){
            imwrite(parameter.debug_dir + "bgMask" + to_string(k) + ".png", U_new_weight_mask[k]);
            U_new_weight_mask[k].copyTo(warp_images_mask[k]);
            /*
            warp_images_mask[k] = false_mask_clean_M(warp_images_mask[k], image_size);
            bitwise_not(warp_images_mask[k],warp_images_mask[k]);
            warp_images_mask[k]=false_mask_clean_M(warp_images_mask[k], image_size);
            bitwise_not(warp_images_mask[k],warp_images_mask[k]);
            imwrite(parameter.debug_dir + "bgMask_optimized" + to_string(k) + ".png", warp_images_mask[k]);
            */
        }
        
    }

    return Blending(_warp_images, warp_images_mask, origins, _target_size, new_weight_mask, _blend_method != BLEND_LINEAR);
}

void MultiImages::writeResultWithMesh(const Mat & _result,
                                      const vector<vector<Point2> > & _vertices,
                                      const string & _postfix,
                                      const bool _only_border) const {
    const int line_thickness = 2;
    const Mat result(_result.size() + Size(line_thickness * 6, line_thickness * 6), CV_8UC4);
    const Point2 shift(line_thickness * 3, line_thickness * 3);
    const Rect rect(shift, _result.size());
    _result.copyTo(result(rect));
    for(int i = 0; i < images_data.size(); ++i) {
        const Scalar & color = getBlueToRedScalar((2. * i / (images_data.size() - 1)) - 1) * 255;
        const vector<Edge> & edges = images_data[i].mesh_2d->getEdges();
        vector<int> edge_indices;
        if(_only_border) {
            edge_indices = images_data[i].mesh_2d->getBoundaryEdgeIndices();
        } else {
            edge_indices.reserve(edges.size());
            for(int j = 0; j < edges.size(); ++j) {
                edge_indices.emplace_back(j);
            }
        }
        for(int j = 0; j < edge_indices.size(); ++j) {
            line(result,
                 _vertices[i][edges[edge_indices[j]].indices[0]] + shift,
                 _vertices[i][edges[edge_indices[j]].indices[1]] + shift, color, line_thickness, LINE_8);
        }
    }
    
    imwrite(parameter.debug_dir + parameter.file_name + _postfix + ".png", result(rect));
}

void MultiImages::writeImageOfFeaturePairs(const string & _name,
                                           const pair<int, int> & _match_pair,
                                           const vector<pair<int, int> > & _pairs) const {
    cout << images_data[_match_pair.first ].file_name << "-" <<
            images_data[_match_pair.second].file_name << " " << _name << " feature pairs = " << _pairs.size() << endl;
    
    const vector<Point2> & m1_fpts = images_data[_match_pair.first ].getFeaturePoints();
    const vector<Point2> & m2_fpts = images_data[_match_pair.second].getFeaturePoints();
    vector<Point2> f1, f2;
    f1.reserve(_pairs.size());
    f2.reserve(_pairs.size());
    for(int i = 0; i < _pairs.size(); ++i) {
        f1.emplace_back(m1_fpts[_pairs[i].first ]);
        f2.emplace_back(m2_fpts[_pairs[i].second]);
    }
    Mat image_of_feauture_pairs = getImageOfFeaturePairs(images_data[_match_pair.first ].img,
                                                         images_data[_match_pair.second].img,
                                                         f1, f2);

    imwrite(parameter.debug_dir +
            "feature_pairs-" + _name + "-" +
            images_data[_match_pair.first ].file_name  + "-" +
            images_data[_match_pair.second].file_name  + "-" +
            to_string(_pairs.size()) +
            images_data[_match_pair.first ].file_extension, image_of_feauture_pairs);
}

vector<vector<vector<InterpolateVertex> > >  MultiImages::get_lsd_mesh_intersection(vector<vector<pair<Point2f, Point2f>>> &endpoints) const{
    vector<vector<vector<InterpolateVertex>>> res;
    res.resize(images_data.size());
    for(int k = 0; k < images_data.size(); k++){
        double lh = images_data[k].mesh_2d->lh;
		double lw = images_data[k].mesh_2d->lw;
		double diag = sqrt(lh*lh + lw*lw);
        vector<double> image_lines_width, image_lines_length,image_lines_slope;
        vector<Vec4f> image_lines;
        vector<Point2f> lines_points[2];
		vector<double> lines_width, lines_length, lines_slope;
        vector<LineData> lines_data = images_data[k].getLines();

        endpoints.resize(images_data.size());
		for (int i = 0; i < lines_data.size(); i++)
		{
		    endpoints[k].push_back(make_pair(lines_data[i].data[0], lines_data[i].data[1]));
		}

        const int lines_count = lines_data.size();
		lines_slope.reserve(lines_count);

        for (int i = 0; i < lines_count; ++i) {
		    lines_width.emplace_back(lines_data[i].width);
		    lines_length.emplace_back(lines_data[i].length);
		    lines_slope.emplace_back((lines_data[i].data[1].y - lines_data[i].data[0].y)/ (lines_data[i].data[1].x - lines_data[i].data[0].x));
		}
        double mean_length = accumulate(begin(lines_length), end(lines_length), 0.0) / lines_length.size();
        double mean_width = accumulate(begin(lines_width), end(lines_width), 0.0) / lines_width.size();
        for(int j = 0; j < lines_count; j++){
            if ((lines_width[j] > 0.5*mean_width) && (lines_length[j] > 0.5*mean_length)&&(lines_length[j]>diag))
                {
                    image_lines.emplace_back(lines_data[j].data[0].x, lines_data[j].data[0].y, lines_data[j].data[1].x, lines_data[j].data[1].y);
                    image_lines_width.push_back(lines_width[j]);
                    image_lines_length.push_back(lines_length[j]);
                    image_lines_slope.push_back(lines_slope[j]);
                }
        }
        res[k].resize(image_lines.size());

        for(int j = 0; j < res[k].size(); j++){
            double x = image_lines[j][0];
            double y = image_lines[j][1];
            int sample_points_num = 0;
            vector <Point2d> sample_points;
            for (int i=0; x <= image_lines[j][2] && x>=image_lines[j][0]; i++){
                sample_points.push_back(Point2d(x, y));
				x += (sqrt(2)*lw*cos(atan(image_lines_slope[j])));
				y += (sqrt(2)*lw*sin(atan(image_lines_slope[j])));
				float min_y = (image_lines[j][3] > image_lines[j][1]) ? image_lines[j][1] : image_lines[j][3];
				float max_y= (image_lines[j][3] > image_lines[j][1]) ? image_lines[j][3] : image_lines[j][1];
				if(y>max_y||y<min_y)
					y -= (sqrt(2)*lw * sin(atan(image_lines_slope[j])));
				sample_points_num ++;    
            }

            res[k][j].reserve(sample_points_num);
            for (int ii = 0; ii < sample_points.size(); ii++)
				res[k][j].emplace_back(images_data[k].mesh_2d->getInterpolateVertex(sample_points[ii]));
        }
    }
    return res;
}

int MultiImages::get_all_line_match_count() const {

    if(all_line_match_count == 0 && line_match_result.empty())
        line_match_result = get_line_match_result();

    return all_line_match_count;
}

vector<vector<vector<pair<pair<cv::Point,cv::Point>,pair<cv::Point,cv::Point>>>>> MultiImages::get_line_match_result() const{
    if(line_match_result.empty()){
        all_line_match_count = 0;
        line_match_result.resize(images_data.size());
        for(int i = 0; i < images_data.size(); i++){
            line_match_result[i].resize(images_data.size());
        }
        vector<vector<vector<KeyLine>>> line_feature_matches = getLineFeatureMatches();
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        for(int i = 0; i < images_match_graph_pair_list.size(); ++i){
            const pair<int, int> &match_pair = images_match_graph_pair_list[i];
            int m1 = match_pair.first;
            int m2 = match_pair.second;
            vector<pair<Point2f,Point2f>> lines_pair_pt1,lines_pair_pt2;
            //if(overlap_rect.empty())
                //overlap_rect = get_overlap_rect();//获取两幅图片重叠区域
            //if(fun.empty())
            //    fun = getFun();//获取两幅图像基础矩阵
            vector<pair<pair<cv::Point,cv::Point>,pair<cv::Point,cv::Point>>> lines_match;
            for(int i=0; i < line_feature_matches[m1][m2].size(); i++){
                lines_pair_pt1.emplace_back(make_pair(line_feature_matches[m1][m2][i].getStartPoint(), line_feature_matches[m1][m2][i].getEndPoint()));
                lines_pair_pt2.emplace_back(make_pair(line_feature_matches[m2][m1][i].getStartPoint(), line_feature_matches[m2][m1][i].getEndPoint()));
            }
            assert(lines_pair_pt1.size() == lines_pair_pt2.size());
            for(int j=0; j < lines_pair_pt1.size(); j++){
                lines_match.emplace_back(make_pair(lines_pair_pt1[j], lines_pair_pt2[j]));
            }
            line_match_result[m1][m2] = lines_match;
            all_line_match_count += lines_match.size();
            //cout<<"all_line_match_count: "<<all_line_match_count<<endl;
            //cout<<"line_match_result[m1][m2].size(): "<<line_match_result[m1][m2].size()<<endl;
        }
    }
    return line_match_result;
}

vector<vector<Rect2d>> MultiImages::get_overlap_rect() const {
    if(overlap_rect.empty())
        getFeaturePairs();
    return overlap_rect;
}

vector<vector<Mat > > MultiImages::getFun() const
{
    if(fun.empty()) {
        if (feature_matches.empty())
            feature_matches = getFeatureMatches();
        const vector<pair<int, int> > &images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        fun.resize(images_data.size());
        for (int i = 0; i < images_data.size(); ++i) {
            fun[i].resize(images_data.size());
        }
        for (int i = 0; i < images_match_graph_pair_list.size(); ++i) {
            const pair<int, int> &match_pair = images_match_graph_pair_list[i];
            const int &m1 = match_pair.first, &m2 = match_pair.second;
            Mat fundamental = findFundamentalMat(feature_matches[m1][m2], feature_matches[m2][m1], CV_FM_RANSAC);
            fun[m1][m2].push_back(fundamental);
            fundamental = findFundamentalMat(feature_matches[m2][m1], feature_matches[m1][m2], CV_FM_RANSAC);
            fun[m2][m1].push_back(fundamental);
        }
    }
    return fun;
}

vector<vector<vector<pair<Point2f,Point2f>>>> MultiImages::get_line_match_sample_points() const{
    if(line_match_sample_points.empty()){
        line_match_sample_points.resize(images_data.size());
        for(int i = 0; i < images_data.size(); i++){
            line_match_sample_points[i].resize(images_data.size());
        }
        const vector<pair<int, int> > & images_match_graph_pair_list = parameter.getImagesMatchGraphPairList();
        vector<vector<vector<pair<pair<cv::Point,cv::Point>,pair<cv::Point,cv::Point>>>>> match_lines = get_line_match_result();
        for(int i = 0; i < images_match_graph_pair_list.size(); i++){
            const pair<int, int> & match_pair = images_match_graph_pair_list[i];
            //cout<< "match_lines[match_pair.first][match_pair.second].size(): "<<match_lines[match_pair.first][match_pair.first].size()<<endl;
            for(int j = 0; j < match_lines[match_pair.first][match_pair.second].size(); j++){
                //pt1、pt5分别为match_pair.first和match_pair.second幅图像间第j对匹配直线在match_pair.first中的两个顶点
                //得到两个顶点后再用二分法均匀地采三个中间点
                Point pt1=match_lines[match_pair.first][match_pair.second][j].first.first;
                Point pt5=match_lines[match_pair.first][match_pair.second][j].first.second;
                Point pt3=0.5*(pt1+pt5);
                Point pt2=0.5*(pt1+pt3);
                Point pt4=0.5*(pt3+pt5);

                //获取这一对匹配直线在match_pair.second中的两个顶点并采样
                Point ptt1=match_lines[match_pair.first][match_pair.second][j].second.first;
                Point ptt5=match_lines[match_pair.first][match_pair.second][j].second.second;
                Point ptt3=0.5*(ptt1+ptt5);
                Point ptt2=0.5*(ptt1+ptt3);
                Point ptt4=0.5*(ptt3+ptt5);
                //将采样得到的五对点保存到line_match_sample_points对象中
                line_match_sample_points[match_pair.first][match_pair.second].emplace_back(make_pair(pt1,ptt1));
                line_match_sample_points[match_pair.first][match_pair.second].emplace_back(make_pair(pt2,ptt2));
                line_match_sample_points[match_pair.first][match_pair.second].emplace_back(make_pair(pt3,ptt3));
                line_match_sample_points[match_pair.first][match_pair.second].emplace_back(make_pair(pt4,ptt4));
                line_match_sample_points[match_pair.first][match_pair.second].emplace_back(make_pair(pt5,ptt5));
            }
        }
    }
    return line_match_sample_points;
}