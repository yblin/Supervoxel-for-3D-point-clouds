//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin.xmu@qq.com (Yangbin Lin)
//

#ifndef VCCS_KNN_SUPERVOXEL_H_
#define VCCS_KNN_SUPERVOXEL_H_

#include <algorithm>
#include <vector>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/macros.h"
#include "codelibrary/geometry/kernel/point_3d.h"
#include "codelibrary/geometry/point_cloud/pca_estimate_normals.h"
#include "codelibrary/geometry/util/distance_3d.h"
#include "codelibrary/util/tree/kd_tree.h"
#include "codelibrary/util/tree/octree.h"

/// KNN variant of VCCS.
class VCCSKNNSupervoxel {
public:
    /**
     * A voxel represents a subset of point cloud, in a regularly spaced,
     * three-dimensional grid.
     */
    struct Voxel {
        int cnt;
        cl::RPoint3D centroid;      /// Centroid point of this voxel.
    };

    /**
     * The supervoxel is represented by a pair (point, normal).
     */
    struct Supervoxel {
        Supervoxel() {}

        Supervoxel(const cl::RPoint3D& p, const cl::RVector3D& n)
            : point(p), normal(n) {}

        cl::RPoint3D point;
        cl::RVector3D normal;
    };

    /**
     * Construct the supervoxel by given [first, last) point cloud and the
     * resolutions.
     */
    VCCSKNNSupervoxel(const cl::KDTree<cl::RPoint3D>& kd_tree,
                      double seed_resolution)
        : seed_resolution_(seed_resolution),
          spatial_importance_(0.4),
          normal_importance_(1.0),
          kd_tree_(kd_tree) {
        assert(seed_resolution_ > 0.0);

        size_points_ = kd_tree_.size();

        InitialSupervoxelSeeds(kd_tree_.points().begin(),
                               kd_tree_.points().end());

        neighbors_.resize(size_points_);
        normals_.resize(size_points_);
        for (int i = 0; i < kd_tree_.size(); ++i) {
            kd_tree_.FindKNearestNeighbors(kd_tree_.points()[i], 20,
                                           &neighbors_[i]);
        }
        cl::Array<cl::RPoint3D> neighbor_points;
        for (int i = 0; i < kd_tree_.size(); ++i) {
            neighbor_points.resize(neighbors_[i].size());
            int k = 0;
            for (int j : neighbors_[i]) {
                neighbor_points[k++] = kd_tree_.points()[j];
            }

            cl::geometry::point_cloud::PCAEstimateNormal(
                        neighbor_points.begin(), neighbor_points.end(),
                        &normals_[i]);
        }
    }

    void set_spatial_importance(double spatial_importance) {
        assert(spatial_importance >= 0.0);

        spatial_importance_ = spatial_importance;
    }

    void set_normal_importance(double normal_importance) {
        assert(normal_importance >= 0.0);

        normal_importance_ = normal_importance;
    }

    /**
     * Segment the given point cloud. It returns an array 'labels' and
     * supervoxels.
     * labels[i] denotes the id of supervoxel that owns the i-th point.
     */
    void Segment(cl::Array<int>* labels,
                 cl::Array<Supervoxel>* supervoxels) {
        assert(labels);
        assert(supervoxels);

        labels->resize(size_points_);
        std::fill(labels->begin(), labels->end(), -1);
        supervoxels->resize(seeds_.size());

        cl::Array<int> queue = seeds_;
        cl::Array<double> distances(size_points_, DBL_MAX);

        for (int i = 0; i < seeds_.size(); ++i) {
            distances[seeds_[i]] = 0.0;
            (*labels)[seeds_[i]] = i;
        }

        const cl::Array<cl::RPoint3D>& points = kd_tree_.points();

        while (!queue.empty()) {
            cl::Array<int> new_queue;

            for (int cur : queue) {
                int label = (*labels)[cur];
                const Supervoxel& facet = (*supervoxels)[(*labels)[cur]];

                for (int neighbor : neighbors_[cur]) {
                    if ((*labels)[neighbor] == label) continue;

                    double dis = FacetDistance(neighbor, facet);
                    if (dis < distances[neighbor]) {
                        distances[neighbor] = dis;
                        (*labels)[neighbor] = label;
                        new_queue.push_back(neighbor);
                    }
                }
            }
            if (new_queue.empty()) break;
            queue = new_queue;

            // Update facets.
            cl::Array<cl::Array<int> > clusters(seeds_.size());
            for (int i = 0; i < size_points_; ++i) {
                int label = (*labels)[i];
                if (label == -1) continue;

                clusters[label].push_back(i);
            }

            for (int i = 0; i < seeds_.size(); ++i) {
                cl::Array<int>& cluster = clusters[i];
                if (cluster.empty()) continue;

                Supervoxel& facet = (*supervoxels)[i];

                cl::RPoint3D centroid;
                cl::RVector3D normal;
                for (int id : cluster) {
                    centroid.x += points[id].x;
                    centroid.y += points[id].y;
                    centroid.z += points[id].z;

                    cl::RVector3D v = normals_[id];
                    if (normal * v > 0.0) {
                        normal += v;
                    } else {
                        normal += -v;
                    }
                }
                normal *= 1.0 / normal.norm();

                centroid.x /= cluster.size();
                centroid.y /= cluster.size();
                centroid.z /= cluster.size();

                double distance = DBL_MAX;
                int nearest_point = -1;
                for (int id : cluster) {
                    double dis = cl::geometry::Distance(points[id], centroid);
                    if (dis < distance) {
                        distance = dis;
                        nearest_point = id;
                    }
                }

                facet = Supervoxel(points[nearest_point], normal);
            }
        }
    }

    void clear() {
        neighbors_.clear();
        neighbors_.shrink_to_fit();
        normals_.clear();
        normals_.shrink_to_fit();
    }

private:
    /**
     * Initial supervoxel seeds.
     */
    template <typename Iterator>
    void InitialSupervoxelSeeds(Iterator first, Iterator last) {
        cl::Array<cl::RPoint3D> seed_points;

        cl::RBox3D box(first, last);
        int size1 = box.x_length() / seed_resolution_ + 1;
        int size2 = box.y_length() / seed_resolution_ + 1;
        int size3 = box.z_length() / seed_resolution_ + 1;

        cl::Octree<int> octree(size1, size2, size3);
        typedef typename cl::Octree<int>::LeafNode LeafNode;

        cl::Array<Voxel> voxels;

        // Add the voxels into the octree.
        int index = 0;
        for (Iterator p = first; p != last; ++p, ++index) {
            int x = (p->x - box.x_min()) / seed_resolution_;
            int y = (p->y - box.y_min()) / seed_resolution_;
            int z = (p->z - box.z_min()) / seed_resolution_;
            x = cl::Clamp(x, 0, size1 - 1);
            y = cl::Clamp(y, 0, size2 - 1);
            z = cl::Clamp(z, 0, size3 - 1);

            std::pair<LeafNode*, bool> pair =
                    octree.Insert(x, y, z, voxels.size());
            if (pair.second) {
                Voxel v;
                v.cnt = 0;
                voxels.push_back(v);
            }
            Voxel& voxel = voxels[pair.first->data()];
            ++voxel.cnt;
            voxel.centroid.x += p->x;
            voxel.centroid.y += p->y;
            voxel.centroid.z += p->z;
        }

        // Compute the adjacent voxels for each voxel in the Octree.
        for (const Voxel& v : voxels) {
            seed_points.emplace_back(v.centroid.x / v.cnt,
                                     v.centroid.y / v.cnt,
                                     v.centroid.z / v.cnt);
        }

        seeds_.resize(seed_points.size());
        for (int i = 0; i < seed_points.size(); ++i) {
            kd_tree_.FindNearestPoint(seed_points[i], &seeds_[i]);
        }
    }

    /**
     * Compute the distance from i-th point to facet.
     */
    double FacetDistance(int i, const Supervoxel& facet) const {
        double t = normals_[i] * facet.normal;
        double n_dist = 1.0 - std::fabs(t);
        double s_dist = cl::geometry::Distance(kd_tree_.points()[i],
                                               facet.point) /
                        seed_resolution_;
        return normal_importance_ * n_dist + spatial_importance_ * s_dist;
    }

    // Number of points in the input point cloud.
    int size_points_;

    // Resolution used to seed the supervoxels.
    double seed_resolution_;

    // Importance of distance from seed center in clustering.
    double spatial_importance_;

    // Importance of similarity in normals for clustering.
    double normal_importance_;

    // The seed supervoxels.
    cl::Array<int> seeds_;

    cl::Array<cl::Array<int> > neighbors_;

    cl::Array<cl::RVector3D> normals_;

    const cl::KDTree<cl::RPoint3D>& kd_tree_;
};

#endif // VCCS_KNN_SUPERVOXEL_H_
