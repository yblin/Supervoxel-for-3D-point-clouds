FROM ubuntu:22.10
RUN apt update
RUN apt install -y build-essential cmake
COPY . code
WORKDIR code
RUN mkdir build
WORKDIR build
RUN cmake .. && cmake --build .
# RUN gcc main.cc -o supervoxel