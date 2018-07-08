//
// Copyright 2013 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_TREE_KD_TREE_H_
#define UTIL_TREE_KD_TREE_H_

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>

#include "codelibrary/base/algorithm.h"
#include "codelibrary/base/array.h"
#include "codelibrary/base/macros.h"
#include "codelibrary/base/object_pool.h"
#include "codelibrary/util/metric/squared_euclidean.h"

namespace cl {

/// KD Tree.
/**
 * KD tree is a space-partitioning data structure for organizing points in a
 * k-dimensional space.
 */
template <typename Point, typename metric = metric::SquaredEuclidean>
class KDTree {
    typedef typename Point::value_type T;

    // The maximum number of points at leaf.
    static const int MAX_LEAF_SIZE = 10;

    /**
     * Bounding box.
     */
    struct BoundingBox {
        BoundingBox() {}

        template <typename Iterator>
        BoundingBox(Iterator first, Iterator last) {
            if (first == last) return;

            int size = first->size();
            assert(size >= 0);

            min_values.resize(size);
            max_values.resize(size);

            Iterator iter = first;
            for (int i = 0; i < size; ++i) {
                min_values[i] = max_values[i] = (*iter)[i];
            }

            for (++iter; iter != last; ++iter) {
                assert(iter->size() == size);

                for (int i = 0; i < size; ++i) {
                    min_values[i] = std::min(min_values[i], (*iter)[i]);
                    max_values[i] = std::max(max_values[i], (*iter)[i]);
                }
            }
        }

        Array<T> min_values, max_values;
    };

    /**
     * A result-set class used when performing a k-nearest Neighbors based
     * search.
     */
    struct KNNResultSet {
        explicit KNNResultSet(int k, double radius = DBL_MAX)
            : capacity(k),
              count(0),
              indices(k, -1),
              distances(k, radius) {}

        /**
         * If the distance from the query point to the given point is smaller
         * than the distances to the other exist points, then add this point
         * into result set.
         *
         * @param index    - the index of point in the KD tree.
         * @param distance - the distance from the query point to this point.
         */
        void AddPoint(int index, double distance) {
            if (distance >= distances.back()) return;

            int i;
            for (i = count; i > 0; --i) {
                if (distance < distances[i - 1]) {
                    if (i < capacity) {
                        distances[i] = distances[i - 1];
                        indices[i] = indices[i - 1];
                    }
                } else {
                    break;
                }
            }

            if (i < capacity) {
                distances[i] = distance;
                indices[i] = index;
            }

            if (count < capacity) ++count;
        }

        /**
         * Return the current farthest distance in the result set.
         */
        double farthest_distance() const {
            return distances.back();
        }

        int capacity;
        int count;
        Array<int> indices;
        Array<double> distances;
    };

    /**
     * A result-set class used when performing a radius based search.
     */
    struct RadiusResultSet {
        explicit RadiusResultSet(double fixed_radius)
            : radius(fixed_radius) {}

        /**
         * If the distance from the query point to the given point is smaller
         * than a fixed radius, then add this point into result set.
         *
         * @param index    - the index of point in the KD tree.
         * @param distance - the distance from the query point to this point.
         */
        void AddPoint(int index, double distance) {
            if (distance <= radius) {
                indices.push_back(index);
            }
        }

        /**
         * Return the farthest distance in result set.
         */
        double farthest_distance() const {
            return radius;
        }

        double radius;
        Array<int> indices;
    };

public:
    /**
     * The Node for KDTree.
     *
     * We public the node structure to help acess the KD tree.
     */
    struct Node {
        union {
            struct {
                int left, right; // Indices of points in leaf node.
            } lr;
            struct {
                int div_dimension;   // Dimension used for subdivision.
                T div_low, div_high; // The values used for subdivision.
            } sub;
        };
        Node* left_child;  // Left child pointer.
        Node* right_child; // Right child pointer.
    };

    /**
     * KDTree constructor.
     */
    KDTree()
        : size_(0), dimension_(0), root_node_(NULL) {}

    template <typename Iterator>
    KDTree(Iterator first, Iterator last)
        : dimension_(0),
          root_node_(NULL),
          points_(first, last) {
        size_ = CountElements(first, last);

        Build();
    }

    explicit KDTree(const Array<Point>& points)
        : KDTree(points.begin(), points.end()) {}

    /**
     * Reset the input points and clear the current KD tree.
     */
    void ResetPoints(const Array<Point>& points) {
        clear();
        points_ = points;
        size_ = points_.size();

        Build();
    }

    /**
     * Swap the KD tree points and rebuild KD tree.
     */
    void SwapPoints(Array<Point>* points) {
        assert(points);

        points_.swap(*points);
        size_ = points_.size();

        Build();
    }

    /**
     * Clear the KD-tree.
     */
    void clear() {
        root_node_ = NULL;
        size_ = 0;
        dimension_ = 0;
        points_.clear();
        indices_.clear();
        node_pool_.clear();
    }

    /**
     * Find the index of nearest neighbor point in the KD tree to the given
     * point.
     */
    void FindNearestPoint(const Point& p, int* nearest_neighbor) const {
        assert(nearest_neighbor);

        // KD tree can not be empty.
        assert(!empty());

        Array<double> distances(dimension_, 0.0);
        double distance_sqr;
        ComputeInitialDistances(p, &distance_sqr, &distances);

        KNNResultSet results(1);
        SearchLevel(root_node_, p, distance_sqr, &distances, &results);
        assert(results.count == 1);

        *nearest_neighbor = results.indices[0];
    }

    /**
     * Find the nearest neighbor point in the KD tree to the given point.
     */
    void FindNearestPoint(const Point& p, Point* nearest_neighbor) const {
        assert(nearest_neighbor);

        int index;
        FindNearestPoint(p, &index);
        *nearest_neighbor = points_[index];
    }

    /**
     * Find the k-nearest neighbors of the given point.
     */
    void FindKNearestNeighbors(const Point& p, int k,
                               Array<int>* neighbors) const {
        assert(neighbors);
        assert(k > 0 && k <= size_);

        Array<double> distances(dimension_, 0.0);
        double distance_sqr;
        ComputeInitialDistances(p, &distance_sqr, &distances);

        KNNResultSet results(k);
        SearchLevel(root_node_, p, distance_sqr, &distances, &results);

        neighbors->swap(results.indices);
        assert(neighbors->size() == k);
    }

    /**
     * Find the k-nearest neighbors of the given point.
     */
    void FindKNearestNeighbors(const Point& p, int k,
                               Array<Point>* neighbors) const {
        assert(neighbors);

        Array<int> indices;
        FindKNearestNeighbors(p, k, &indices);

        neighbors->resize(indices.size());
        for (int i = 0; i < indices.size(); ++i) {
            (*neighbors)[i] = points_[indices[i]];
        }
    }

    /**
     * Find the k-nearest neighbors within given radius of the given point.
     */
    void FindKNearestInRadiusNeighbors(const Point& p, int k, double radius,
                                       Array<int>* neighbors) const {
        assert(neighbors);
        assert(k > 0 && k <= size_);
        assert(radius > 0.0);

        Array<double> distances(dimension_, 0.0);
        double distance_sqr;
        ComputeInitialDistances(p, &distance_sqr, &distances);

        KNNResultSet results(k, radius);
        SearchLevel(root_node_, p, distance_sqr, &distances, &results);

        neighbors->swap(results.indices);
        neighbors->resize(results.count);
    }

    /**
     * Find the k-nearest neighbors within given radius of given point.
     */
    void FindKNearestInRadiusNeighbors(const Point& p, int k, double radius,
                                       Array<Point>* neighbors) const {
        assert(neighbors);

        Array<int> indices;
        FindKNearestInRadiusNeighbors(p, k, radius, &indices);

        neighbors->resize(indices.size());
        for (int i = 0; i < indices.size(); ++i) {
            (*neighbors)[i] = points_[indices[i]];
        }
    }

    /**
     * Find the fixed-radius neighbors of given point.
     */
    void FindRadiusNeighbors(const Point& p, double radius,
                             Array<int>* neighbors) const {
        assert(neighbors);
        assert(!empty());

        Array<double> distances(dimension_, 0.0);
        double distance_sqr;
        ComputeInitialDistances(p, &distance_sqr, &distances);

        RadiusResultSet results(radius);
        SearchLevel(root_node_, p, distance_sqr, &distances, &results);

        neighbors->swap(results.indices);
    }

    /**
     * Find the fixed-radius neighbors of the given point.
     */
    void FindRadiusNeighbors(const Point& p, double radius,
                             Array<Point>* neighbors) const {
        assert(neighbors);

        Array<int> indices;
        FindRadiusNeighbors(p, radius, &indices);

        neighbors->resize(indices.size());
        for (int i = 0; i < indices.size(); ++i) {
            (*neighbors)[i] = points_[indices[i]];
        }
    }

    int size()                   const { return size_;      }
    bool empty()                 const { return size_ == 0; }
    const metric& distance()     const { return distance_;  }
    const Array<Point>& points() const { return points_;    }
    const Node* root_node()      const { return root_node_; }

private:
    /**
     * Build the KD tree.
     */
    void Build() {
        if (empty()) return;

        dimension_ = points_[0].size();
        assert(dimension_ > 0);

        node_pool_.clear();
        indices_.resize(size_);
        for (int i = 0; i < size_; ++i) {
            indices_[i] = i;
        }
        box_ = BoundingBox(points_.begin(), points_.end());

        root_node_ = DivideTree(0, size_, box_);
    }

    /**
     * @param[in]  p            - the point for searching.
     * @param[out] distance_sqr - the square of distance between point to
     *                            bounding box.
     * @param[out] distances    - the distance between point to each dimension
     *                            of bounding box (need initialized before).
     */
    void ComputeInitialDistances(const Point& p, double* distance_sqr,
                                 Array<double>* distances) const {
        *distance_sqr = 0.0;

        for (int i = 0; i < dimension_; ++i) {
            if (p[i] < box_.min_values[i]) {
                (*distances)[i] = (p[i] - box_.min_values[i]) *
                                  (p[i] - box_.min_values[i]);
                *distance_sqr += (*distances)[i];
            }
            if (p[i] > box_.max_values[i]) {
                (*distances)[i] = (p[i] - box_.max_values[i]) *
                                  (p[i] - box_.max_values[i]);
                *distance_sqr += (*distances)[i];
            }
        }
    }

    /**
     * Perform an exact search in the tree starting from a node.
     */
    template <typename ResultSet>
    void SearchLevel(const Node* node, const Point& query_point,
                     double min_distance_sqr, Array<double>* distances,
                     ResultSet* results) const {
        // If this this is a leaf node, then do check and return.
        if (node->left_child == NULL && node->right_child == NULL) {
            for (int i = node->lr.left; i < node->lr.right; ++i) {
                double dis = distance_(points_[indices_[i]], query_point);
                if (dis < results->farthest_distance())
                    results->AddPoint(indices_[i], dis);
            }
            return;
        }

        // Find which child branch should be taken first.
        int d = node->sub.div_dimension;
        T value = query_point[d];
        double diff1 = value - node->sub.div_low;
        double diff2 = value - node->sub.div_high;

        Node* best_child;
        Node* other_child;
        double cut_distance;
        if ((diff1 + diff2) < 0.0) {
            best_child = node->left_child;
            other_child = node->right_child;
            double t = value - node->sub.div_high;
            cut_distance = t * t;
        } else {
            best_child = node->right_child;
            other_child = node->left_child;
            double t = value - node->sub.div_low;
            cut_distance = t * t;
        }

        // Call recursively to search next level down.
        SearchLevel(best_child, query_point, min_distance_sqr, distances,
                    results);

        double distance = (*distances)[d];

        // Compute the square of distance between search point to the bounding
        // box of the points in other_child.
        double lower_bound = min_distance_sqr + cut_distance - distance;

        (*distances)[d] = cut_distance;
        if (lower_bound <= results->farthest_distance()) {
            SearchLevel(other_child, query_point, lower_bound, distances,
                        results);
        }
        (*distances)[d] = distance;
    }

    /**
     * Create a tree node that subdivides the list of vertices from
     * indices[left, right).
     *
     * The routine is called recursively on each sublist.
     */
    Node* DivideTree(int left, int right, const BoundingBox& box) {
        Node* node = node_pool_.Allocate();

        // If too few exemplars remain, then make this a leaf node.
        if ((right - left) <= MAX_LEAF_SIZE) {
            node->left_child = node->right_child = NULL; // Mark as leaf node.
            node->lr.left = left;
            node->lr.right = right;
        } else {
            int cut_index;
            int cut_dimension;
            T cut_value;
            MiddleSplit(left, right, box, &cut_index, &cut_dimension,
                        &cut_value);

            node->sub.div_dimension = cut_dimension;

            BoundingBox left_box(box);
            left_box.max_values[cut_dimension] = cut_value;
            node->left_child = DivideTree(left, cut_index, left_box);

            BoundingBox right_box(box);
            right_box.min_values[cut_dimension] = cut_value;
            node->right_child = DivideTree(cut_index, right, right_box);

            node->sub.div_low = left_box.max_values[cut_dimension];
            node->sub.div_high = right_box.min_values[cut_dimension];
        }

        return node;
    }

    /**
     * Compute the minimum and maximum element value of k-th component of
     * points in [left, right).
     */
    void ComputeMinMax(int left, int right, int k, T* min_elem, T* max_elem) {
        *min_elem = *max_elem = points_[indices_[left]][k];
        for (int i = left + 1; i < right; ++i) {
            T val = points_[indices_[i]][k];
            *min_elem = std::min(*min_elem, val);
            *max_elem = std::max(*max_elem, val);
        }
    }

    /**
     * Middle split the vertices from indices[left] to indices[right].
     *
     * @param[in]  left is the first index to split.
     * @param[in]  right is the last index to split.
     * @param[in]  box is the approximate bounding box of points[left, right).
     * @param[out] index in [left, right) is the index of split.
     * @param[out] cut_feat is the dimension to split.
     * @param[out] cut_value is position of split.
     */
    void MiddleSplit(int left, int right, const BoundingBox& box,
                     int* cut_index, int* cut_dimension, T* cut_value) {
        // Find the largest span from the approximate bounding box.
        T max_span = box.max_values[0] - box.min_values[0];
        *cut_dimension = 0;
        for (int i = 1; i < dimension_; ++i) {
            T span = box.max_values[i] - box.min_values[i];
            if (span > max_span) {
                max_span = span;
                *cut_dimension = i;
            }
        }

        // Compute exact span on the found dimension.
        T min_elem, max_elem;
        ComputeMinMax(left, right, *cut_dimension, &min_elem, &max_elem);
        *cut_value = (min_elem + max_elem) / 2;
        max_span = max_elem - min_elem;

        // Check if a dimension of a largest span exists.
        for (int i = 0; i < dimension_; ++i) {
            if (i == *cut_dimension) continue;
            T span = box.max_values[i] - box.min_values[i];
            if (span > max_span) {
                ComputeMinMax(left, right, i, &min_elem, &max_elem);
                span = max_elem - min_elem;
                if (span > max_span) {
                    max_span = span;
                    *cut_dimension = i;
                    *cut_value = (min_elem + max_elem) / 2;
                }
            }
        }

        int lim1, lim2;
        PlaneSplit(left, right, *cut_dimension, *cut_value, &lim1, &lim2);

        int count = right - left;
        if (lim1 > left + count / 2) {
            *cut_index = lim1;
        } else if (lim2 < left + count / 2) {
            *cut_index = lim2;
        } else {
            // If either list is empty, it means that all remaining features are
            // identical. Split in the middle to maintain a balanced tree.
            *cut_index = left + count / 2;
        }
    }

    /**
     * Subdivide the list of points by a plane perpendicular on axe
     * corresponding to the 'cut_dimension' dimension at 'cut_value' position.
     *
     * Out parameters:
     *   points[indices[left .. left+lim1-1]][cut_dimension]     < cut_value;
     *   points[indices[left+lim1 .. leftlim2-1]][cut_dimension] = cut_value;
     *   points[indices[left+lim2 .. right]][cut_dimension]      > cut_value.
     */
    void PlaneSplit(int left, int right, int cut_dimension, const T& cut_value,
                    int* lim1, int* lim2) {
        int l = left, r = right - 1;

        // Move vector indices for left subtree to front of list.
        for (;;) {
            while (l <= r &&
                   points_[indices_[l]][cut_dimension] < cut_value) ++l;
            while (r && l <= r &&
                   points_[indices_[r]][cut_dimension] >= cut_value) --r;
            if (l > r || !r) break;
            std::swap(indices_[l], indices_[r]);
            ++l;
            --r;
        }

        *lim1 = l;
        r = right - 1;
        for (;;) {
            while (l <= r &&
                   points_[indices_[l]][cut_dimension] <= cut_value) ++l;
            while (r && l <= r &&
                   points_[indices_[r]][cut_dimension] > cut_value) --r;
            if (l > r || !r) break;
            std::swap(indices_[l], indices_[r]);
            ++l;
            --r;
        }
        *lim2 = l;
    }

    int size_;                   // The number of points.
    int dimension_;              // The dimension of the KD tree.
    BoundingBox box_;            // The bounding box of points.
    metric distance_;            // The distance metric.
    Node* root_node_;            // The pointer to the root of tree.
    Array<Point> points_;        // The points data of KD tree.
    Array<int> indices_;         // The indices of the points.
    ObjectPool<Node> node_pool_; // The memory pool for nodes of KD tree.

    DISALLOW_COPY_AND_ASSIGN(KDTree);
};

} // namespace cl

#endif // UTIL_TREE_KD_TREE_H_
