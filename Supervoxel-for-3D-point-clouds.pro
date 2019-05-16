TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cc

HEADERS += \
    codelibrary/geometry/point_cloud/grid_sample.h \
    codelibrary/geometry/point_cloud/pca_estimate_normals.h \
    codelibrary/base/algorithm.h \
    codelibrary/base/array.h \
    codelibrary/base/bits.h \
    codelibrary/base/equal.h \
    codelibrary/base/log.h \
    codelibrary/base/macros.h \
    codelibrary/base/message.h \
    codelibrary/base/object_pool.h \
    codelibrary/base/string_printf.h \
    codelibrary/geometry/io/xyz_io.h \
    codelibrary/geometry/kernel/box_2d.h \
    codelibrary/geometry/kernel/box_3d.h \
    codelibrary/geometry/kernel/circle_2d.h \
    codelibrary/geometry/kernel/line_2d.h \
    codelibrary/geometry/kernel/line_3d.h \
    codelibrary/geometry/kernel/plane_3d.h \
    codelibrary/geometry/kernel/point_2d.h \
    codelibrary/geometry/kernel/point_3d.h \
    codelibrary/geometry/kernel/point_nd.h \
    codelibrary/geometry/kernel/polygon_2d.h \
    codelibrary/geometry/kernel/segment_2d.h \
    codelibrary/geometry/kernel/segment_3d.h \
    codelibrary/geometry/kernel/sphere_3d.h \
    codelibrary/geometry/kernel/triangle_2d.h \
    codelibrary/geometry/kernel/triangle_3d.h \
    codelibrary/geometry/transform/project_3d.h \
    codelibrary/geometry/util/center_2d.h \
    codelibrary/geometry/util/center_3d.h \
    codelibrary/geometry/util/compare_2d.h \
    codelibrary/geometry/util/compare_3d.h \
    codelibrary/geometry/util/distance_3d.h \
    codelibrary/geometry/point_octree.h \
    codelibrary/math/matrix/decompose/tridiagonal_decompose.h \
    codelibrary/math/matrix/matrix.h \
    codelibrary/math/matrix/matrix_util.h \
    codelibrary/math/matrix/symmetric_eigen.h \
    codelibrary/math/matrix/triplet_matrix.h \
    codelibrary/math/vector/vector.h \
    codelibrary/math/vector/vector_2d.h \
    codelibrary/math/vector/vector_3d.h \
    codelibrary/math/angle.h \
    codelibrary/math/matrix.h \
    codelibrary/math/vector.h \
    codelibrary/statistics/kernel/covariance.h \
    codelibrary/statistics/kernel/deviation.h \
    codelibrary/statistics/kernel/mean.h \
    codelibrary/statistics/kernel/median.h \
    codelibrary/statistics/kernel/ranking.h \
    codelibrary/statistics/regression/linear_least_squares_fitting.h \
    codelibrary/statistics/principal_component_analysis_2d.h \
    codelibrary/statistics/principal_component_analysis_3d.h \
    codelibrary/util/array/array_2d.h \
    codelibrary/util/metric/angular.h \
    codelibrary/util/metric/cosine.h \
    codelibrary/util/metric/euclidean.h \
    codelibrary/util/metric/jaccard.h \
    codelibrary/util/metric/manhattan.h \
    codelibrary/util/metric/squared_euclidean.h \
    codelibrary/util/set/disjoint_set.h \
    codelibrary/util/tree/kd_tree.h \
    codelibrary/util/tree/octree.h \
    codelibrary/visualization/color/rgb32_color.h \
    codelibrary/geometry/point_cloud/supervoxel_segmentation.h \
    vccs_knn_supervoxel.h \
    vccs_supervoxel.h
