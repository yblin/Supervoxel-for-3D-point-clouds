# Supervoxel for 3D point clouds

## Introduction
We present a simple but effective supervoxel segmentation method for point clouds, which formalizes supervoxel segmentation as a subset selection problem. We develop an heuristic algorithm that utilizes local information to efficiently solve the subset selection problem. The proposed method can produce supervoxels with adaptive resolutions, and dose not rely the selection of seed points. The method is fully tested on three publicly available point cloud segmentation benchmarks, which cover the major point cloud types. The experimental results show that compared with the state-of-the-art supervoxel segmentation methods, the supervoxels extracted using our method preserve the object boundaries and small structures more effectively, which is reflected in a higher boundary recall and lower under-segmentation error.

The details can be found in the following [ISPRS 2018 paper](https://www.sciencedirect.com/science/article/pii/S0924271618301370)

## Citing our work
If you find our works useful in your research, please consider citing:

   Lin Y, Wang C, Zhai D, et al. Toward better boundary preserved supervoxel segmentation for 3D point clouds. Isprs Journal of Photogrammetry & Remote Sensing, 2018.
   
## Install & complie

Please directly copy the code into your workspace and complie it with any complier that supports C++11. It dose not require linking any additional libraries.

## Sample usage:
	cl::geometry::point_cloud::SupervoxelSegmentation(oriented_points,
													  neighbors,
													  resolution,
												      metric,
													  &supervoxels,
													  &labels);

Where, 'oriented_points' is the source points, it can be read from XYZ file by calling: 
cl::geometry::io::ReadXYZPoints(filename.c_str(), &points);

See main.cc for more details.

