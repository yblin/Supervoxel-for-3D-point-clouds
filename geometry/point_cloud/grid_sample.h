//
// Copyright 2016 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef GEOMETRY_POINT_CLOUD_GRID_SAMPLE_H_
#define GEOMETRY_POINT_CLOUD_GRID_SAMPLE_H_

#include <algorithm>
#include <limits>
#include <random>

#include "codelibrary/base/array.h"
#include "codelibrary/base/algorithm.h"
#include "codelibrary/geometry/kernel/box_3d.h"
#include "codelibrary/util/tree/octree.h"

namespace cl {
namespace geometry {
namespace point_cloud {

/**
 * It considers a regular grid covering the bounding box of the input point
 * cloud, and clusters all points sharing the same cell of the grid by picking
 * as representant one arbitrarily chosen point.
 */
template <typename Iterator>
void GridSample(Iterator first, Iterator last, double resolution,
                Array<int>* sampling) {
    assert(resolution > 0.0);
    assert(sampling);

    sampling->clear();

    int n = CountElements(first, last);
    assert(n > 0);

    Array<int> random_seq(n);
    for (int i = 0; i < n; ++i) {
        random_seq[i] = i;
    }
    std::mt19937 random;
    std::shuffle(random_seq.begin(), random_seq.end(), random);

    RBox3D box(first, last);
    assert(box.x_length() / resolution < INT_MAX);
    assert(box.y_length() / resolution < INT_MAX);
    assert(box.z_length() / resolution < INT_MAX);
    int size1 = box.x_length() / resolution + 1;
    int size2 = box.y_length() / resolution + 1;
    int size3 = box.z_length() / resolution + 1;

    Octree<bool> octree(size1, size2, size3);
    typedef typename Octree<bool>::LeafNode LeafNode;
    typedef typename std::iterator_traits<Iterator>::value_type Point;

    // Add the voxels into the octree.
    for (int s : random_seq) {
        const Point& p = first[s];
        int x = (p.x - box.x_min()) / resolution;
        int y = (p.y - box.y_min()) / resolution;
        int z = (p.z - box.z_min()) / resolution;
        x = Clamp(x, 0, size1 - 1);
        y = Clamp(y, 0, size2 - 1);
        z = Clamp(z, 0, size3 - 1);

        std::pair<LeafNode*, bool> pair = octree.Insert(x, y, z, true);
        if (pair.second) {
            sampling->push_back(s);
        }
    }
}

} // namespace point_cloud
} // namespace geometry
} // namespace cl

#endif // GEOMETRY_POINT_CLOUD_GRID_SAMPLE_H_
