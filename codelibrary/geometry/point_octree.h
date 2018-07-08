//
// Copyright 2018 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_POINT_OCTREE_H_
#define GEOMETRY_POINT_OCTREE_H_

#include <cassert>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/util/metric/squared_euclidean.h"
#include "codelibrary/util/tree/octree.h"

namespace cl {
namespace geometry {

/**
 * Octree for 3D point set.
 */
template <typename Point, typename Metric = metric::SquaredEuclidean>
class PointOctree {
    typedef typename Octree<Array<int> >::LeafNode LeafNode;

public:
    PointOctree(double resolution)
        : size_(0), resolution_(resolution) {
        assert(resolution_ > 0.0);
    }

    template <typename Iterator>
    PointOctree(Iterator first, Iterator last, double resolution)
        : points_(first, last), resolution_(resolution) {
        assert(resolution_ > 0.0);

        size_ = CountElements(first, last);
        Build();
    }

    PointOctree(const Array<Point>& points, double resolution)
        : PointOctree(points.begin(), points.end()) {}

    /**
     * Reset the point set, and rebuild the octree.
     */
    void ResetPoints(const Array<Point>& points) {
        points_ = points;
        size_ = points_.size();

        Build();
    }

    bool empty() const {
        return size_ == 0;
    }

    /**
     * Find the fixed-radius neighbors of the given point.
     */
    void FindNeighbors(const Point& p, double radius,
                       Array<int>* neighbors) const {
        assert(neighbors);
        assert(!empty());

        neighbors->clear();

        int x = (p.x - bounding_box_.x_min()) / resolution_;
        int y = (p.y - bounding_box_.y_min()) / resolution_;
        int z = (p.z - bounding_box_.z_min()) / resolution_;
        x = Clamp(x, 0, octree_.size1() - 1);
        y = Clamp(y, 0, octree_.size2() - 1);
        z = Clamp(z, 0, octree_.size3() - 1);

        for (int x1 = x - 1; x1 <= x + 1; ++x1) {
            if (x1 < 0 || x1 >= octree_.size1()) continue;
            for (int y1 = y - 1; y1 <= y + 1; ++y1) {
                if (y1 < 0 || y1 >= octree_.size2()) continue;
                for (int z1 = z - 1; z1 <= z + 1; ++z1) {
                    if (z1 < 0 || z1 >= octree_.size3()) continue;

                    const LeafNode* leaf = octree_.Find(x1, y1, z1);
                    if (!leaf) continue;

                    for (int i : leaf->data()) {
                        double dis = distance_(points_[i], p);
                        if (dis <= radius) {
                            neighbors->push_back(points_[i]);
                        }
                    }
                }
            }
        }
    }

    /**
     * Find the fixed-radius neighbors of given point.
     */
    void FindNeighbors(const Point& p, double radius,
                       Array<Point>* neighbors) const {
        assert(neighbors);

        Array<int> indices;
        FindRadiusNeighbors(p, radius, &indices);

        neighbors->resize(indices.size());
        for (int i = 0; i < indices.size(); ++i) {
            (*neighbors)[i] = points_[indices[i]];
        }
    }

    /**
     * Swap the points and rebuild octree.
     */
    void SwapPoints(Array<Point>* points) {
        assert(points);

        points_.swap(*points);
        size_ = points_.size();

        Build();
    }

    /**
     * @return the number of points.
     */
    int size() const {
        return size_;
    }

    /**
     * @return the resolution of octree.
     */
    double resolution() const {
        return resolution_;
    }

    const Metric& distance() const {
        return distance_;
    }

    const Array<Point>& points() const {
        return points_;
    }

private:
    /**
     * Build octree.
     */
    void Build() {
        bounding_box_ = RBox3D(points_.begin(), points_.end());
        assert(bounding_box_.x_length() / resolution_ < INT_MAX);
        assert(bounding_box_.y_length() / resolution_ < INT_MAX);
        assert(bounding_box_.z_length() / resolution_ < INT_MAX);

        int size1 = bounding_box_.x_length() / resolution_ + 1;
        int size2 = bounding_box_.y_length() / resolution_ + 1;
        int size3 = bounding_box_.z_length() / resolution_ + 1;

        octree_.Resize(size1, size2, size3);

        // Add the voxels into the octree.
        Array<int> tmp;
        for (int i = 0; i < size_; ++i) {
            const Point& p = points_[i];
            int x = (p.x - bounding_box_.x_min()) / resolution_;
            int y = (p.y - bounding_box_.y_min()) / resolution_;
            int z = (p.z - bounding_box_.z_min()) / resolution_;
            x = Clamp(x, 0, size1 - 1);
            y = Clamp(y, 0, size2 - 1);
            z = Clamp(z, 0, size3 - 1);

            std::pair<LeafNode*, bool> pair = octree_.Insert(x, y, z, tmp);
            pair.first->data().push_back(i);
        }
    }

    int size_;                    // The number of points.
    double resolution_;           // Resolution to construct octree.
    RBox3D bounding_box_;         // The bounding box of 3D point set.
    Octree<Array<int> > octree_;  // Octree to store the index of point set.
    Array<Point> points_;         // The input points.
    Metric distance_;             // The distance metric.
};

} // namespace geometry
} // namespace cl

#endif // GEOMETRY_POINT_OCTREE_H_
