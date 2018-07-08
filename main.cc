#include "codelibrary/base/log.h"
#include "codelibrary/geometry/io/xyz_io.h"
#include "codelibrary/geometry/point_cloud/pca_estimate_normals.h"
#include "codelibrary/geometry/point_cloud/supervoxel_segmentation.h"
#include "codelibrary/geometry/util/distance_3d.h"
#include "codelibrary/util/tree/kd_tree.h"

/// Point with Normal.
struct PointWithNormal : cl::RPoint3D {
    PointWithNormal() {}

    cl::RVector3D normal;
};

/**
 * Metric used in VCCS supervoxel segmentation.
 *
 * Reference:
 *   Rusu, R.B., Cousins, S., 2011. 3d is here: Point cloud library (pcl),
 *   IEEE International Conference on Robotics and Automation, pp. 1–4.
 */
class VCCSMetric {
public:
    explicit VCCSMetric(double resolution)
        : resolution_(resolution) {}

    double operator() (const PointWithNormal& p1,
                       const PointWithNormal& p2) const {
        return 1.0 - std::fabs(p1.normal * p2.normal) +
               cl::geometry::Distance(p1, p2) / resolution_ * 0.4;
    }

private:
    double resolution_;
};

int main() {
    LOG_ON(INFO);

    const std::string filename = "test.xyz";

    cl::Array<cl::RPoint3D> points;
    cl::Array<cl::RGB32Color> colors;

    LOG(INFO) << "Reading points from " << filename << "...";
    if (!cl::geometry::io::ReadXYZPoints(filename.c_str(), &points, &colors)) {
        return 0;
    }
    int n_points = points.size();
    LOG(INFO) << n_points << " points are imported.";

    LOG(INFO) << "Building KD tree...";
    cl::KDTree<cl::RPoint3D> kdtree;
    kdtree.SwapPoints(&points);

    const int k_neighbors = 15;
    assert(k_neighbors < n_points);

    LOG(INFO) << "Compute the k-nearest neighbors for each point, and "
                 "estimate the normal vectors...";

    cl::Array<cl::RVector3D> normals(n_points);
    cl::Array<cl::Array<int> > neighbors(n_points);
    cl::Array<cl::RPoint3D> neighbor_points(k_neighbors);
    for (int i = 0; i < n_points; ++i) {
        kdtree.FindKNearestNeighbors(kdtree.points()[i], k_neighbors,
                                     &neighbors[i]);
        for (int k = 0; k < k_neighbors; ++k) {
            neighbor_points[k] = kdtree.points()[neighbors[i][k]];
        }
        cl::geometry::point_cloud::PCAEstimateNormal(neighbor_points.begin(),
                                                     neighbor_points.end(),
                                                     &normals[i]);
    }
    kdtree.SwapPoints(&points);

    LOG(INFO) << "Start supervoxel segmentation...";

    cl::Array<PointWithNormal> oriented_points(n_points);
    for (int i = 0; i < n_points; ++i) {
        oriented_points[i].x = points[i].x;
        oriented_points[i].y = points[i].y;
        oriented_points[i].z = points[i].z;
        oriented_points[i].normal = normals[i];
    }

    const double resolution = 0.5;
    VCCSMetric metric(resolution);
    cl::Array<int> labels, supervoxels;
    cl::geometry::point_cloud::SupervoxelSegmentation(oriented_points,
                                                      neighbors,
                                                      resolution,
                                                      metric,
                                                      &supervoxels,
                                                      &labels);

    int n_supervoxels = supervoxels.size();
    LOG(INFO) << n_supervoxels << " supervoxels computed.";

    std::mt19937 random;
    cl::Array<cl::RGB32Color> supervoxel_colors(n_supervoxels);
    for (int i = 0; i < n_supervoxels; ++i) {
        supervoxel_colors[i] = cl::RGB32Color(random());
    }
    for (int i = 0; i < n_points; ++i) {
        colors[i] = supervoxel_colors[labels[i]];
    }

    if (cl::geometry::io::WriteXYZPoints("out.xyz", points, colors)) {
        LOG(INFO) << "The points are written into out.xyz";
    }

    system("out.xyz");

    return 0;
}
