//
// Copyright 2014 Yangbin Lin. All Rights Reserved.
//
// Author: yblin@jmu.edu.cn (Yangbin Lin)
//
// This file is part of the Code Library.
//

#ifndef UTIL_TREE_OCTREE_H_
#define UTIL_TREE_OCTREE_H_

#include <algorithm>
#include <cassert>
#include <cstring>
#include <utility>

#include "codelibrary/base/array.h"
#include "codelibrary/base/bits.h"
#include "codelibrary/base/macros.h"
#include "codelibrary/base/object_pool.h"

namespace cl {

/// Octree.
/**
 * The class Octree can be used as a 3D array, but with higher space efficiency.
 *
 * Sample usage:
 *
 * Octree<int> o(4, 4, 4);
 * o(1, 2, 3) = 3;
 */
template <typename T>
class Octree {
public:
    /// Enum of node type within the Octree.
    enum NodeType {
        BRANCH_NODE, LEAF_NODE
    };

    /// Basic node of Octree for inherit.
    struct Node {
        Node() {}

        virtual ~Node() {}

        virtual NodeType get_type() const = 0;
    };

    /// Branch Node for Octree.
    class BranchNode : Node {
        friend class Octree;

    public:
        BranchNode() {
            memset(child_nodes_, 0, sizeof(child_nodes_));
        }

        virtual NodeType get_type() const {
            return BRANCH_NODE;
        }

    protected:
        Node* child_nodes_[8]; // The pointers to the children.
    };

    /// Leaf Node for Octree.
    class LeafNode : Node {
        friend class Octree;

    public:
        LeafNode() {}

        virtual NodeType get_type() const {
            return LEAF_NODE;
        }

        const T& data() const {
            return data_;
        }

        T& data() {
            return data_;
        }

        void set_data(const T& data) {
            data_ = data;
        }

    protected:
        T data_; // The data stored in this leaf node.
    };

    Octree()
        : size1_(0), size2_(0), size3_(0) {
        Initialize();
    }

    Octree(int size1, int size2, int size3)
        : size1_(size1), size2_(size2), size3_(size3) {
        assert(size1_ >= 0 && size2_ >= 0 && size3_ >= 0);

        Initialize();
     }

    /**
     * Clear the octree.
     */
    void clear() {
        branch_pool_.clear();
        leaf_pool_.clear();

        root_ = branch_pool_.Allocate();
        size_ = 0;
        depth_ = 0;
        size1_ = 0;
        size2_ = 0;
        size3_ = 0;
    }

    /**
     * @return the size of elements.
     */
    int size() const {
        return size_;
    }

    /**
     * @return the size of the first dimension.
     */
    int size1() const {
        return size1_;
    }

    /**
     * @return the size of the second dimension.
     */
    int size2() const {
        return size2_;
    }

    /**
     * @return the size of third dimension.
     */
    int size3() const {
        return size3_;
    }

    /**
     * @return the root of octree.
     */
    BranchNode* root() const {
        return root_;
    }

    /**
     * @return the depth of octree.
     */
    int depth() const {
        return depth_;
    }

    /**
     * Check if the octree is empty.
     */
    bool empty() const {
        return depth_ == 0;
    }

    /**
     * Resize the octree.
     *
     * @note that, this method will clear the data.
     */
    void Resize(int size1, int size2, int size3) {
        assert(size1 >= 0 && size2 >= 0 && size3 >= 0);

        clear();

        size1_ = size1;
        size2_ = size2;
        size3_ = size3;

        SetDepth();
    }

    /**
     * Find the leaf node at (x, y, z) in octree.
     */
    LeafNode* Find(int x, int y, int z) {
        assert(CheckValid(x, y, z));

        return Find(x, y, z, DepthMask(), root_);
    }

    /**
     * Find the leaf node at (x, y, z) in octree.
     */
    const LeafNode* Find(int x, int y, int z) const {
        assert(CheckValid(x, y, z));

        return Find(x, y, z, DepthMask(), root_);
    }

    /**
     * Insert a leaf node at (x, y, z) in octree. If the leaf node already
     * exists, its data will not be changed.
     *
     * The function returns a pair, with its member pair::first set to an
     * iterator pointing to either the newly inserted element or to the
     * equivalent element already in the set. The pair::second element in the
     * pair is set to true if a new element was inserted or false if an
     * equivalent element already existed.
     */
    std::pair<LeafNode*, bool> Insert(int x, int y, int z, const T& v) {
        assert(CheckValid(x, y, z));

        return Insert(x, y, z, DepthMask(), root_, v);
    }

    /**
     * Erase a leaf node, return false if the node does not exist.
     */
    bool Erase(int x, int y, int z) {
        assert(CheckValid(x, y, z));

        return Erase(x, y, z, DepthMask(), root_);
    }

    T& operator() (int x, int y, int z) {
        std::pair<LeafNode*, bool> leaf = Insert(x, y, z, T());
        return leaf.first->data_;
    }

    const T& operator() (int x, int y, int z) const {
        const LeafNode* leaf = Find(x, y, z);

        return (leaf == NULL) ? T() : leaf->data_;
    }

protected:
    /**
     * Check if (x, y, z) is a valid octree index.
     */
    bool CheckValid(int x, int y, int z) const {
        return x >= 0 && x < size1_ &&
               y >= 0 && y < size2_ &&
               z >= 0 && z < size3_;
    }

    /**
     * Insert a leaf node.
     *
     * Note that, if the leaf node already exists, its data will not be changed,
     */
    std::pair<LeafNode*, bool> Insert(int x, int y, int z, int depth_mask,
                                      BranchNode* branch, const T& data) {
        int child = GetChildIndex(x, y, z, depth_mask);
        Node* node = branch->child_nodes_[child];

        if (node == NULL) {
            if (depth_mask > 1) {
                // If required branch does not exist, create it.
                BranchNode* child_branch = branch_pool_.Allocate();
                branch->child_nodes_[child] =
                        reinterpret_cast<Node*>(child_branch);

                return Insert(x, y, z, depth_mask >> 1, child_branch, data);
            } else {
                // If leaf node at child_index does exist.
                LeafNode* leaf = leaf_pool_.Allocate();
                leaf->data_ = data;
                branch->child_nodes_[child] = reinterpret_cast<Node*>(leaf);
                ++size_;
                return std::make_pair(leaf, true);
            }
        } else {
            // Node exists already.
            switch (node->get_type()) {
            case BRANCH_NODE:
                return Insert(x, y, z, depth_mask >> 1,
                              reinterpret_cast<BranchNode*>(node), data);
            case LEAF_NODE:
                // Found.
                return std::make_pair(reinterpret_cast<LeafNode*>(node), false);
            }
        }

        // Code never go here.
        assert(false && "Unreachable code.");
        return std::make_pair(reinterpret_cast<LeafNode*>(NULL), false);
    }

    /**
     * Find a leaf node and return it.
     * If the leaf already exist, directly return it.
     */
    LeafNode* Find(int x, int y, int z, int depth_mask, BranchNode* branch) {
        int child = GetChildIndex(x, y, z, depth_mask);
        Node* node = branch->child_nodes_[child];
        if (node == NULL) return NULL;

        if (node->get_type() == BRANCH_NODE) {
            return Find(x, y, z, depth_mask >> 1,
                        reinterpret_cast<BranchNode*>(node));
        } else {
            return reinterpret_cast<LeafNode*>(node);
        }
    }

    /**
     * Find a leaf node and return it.
     * If the leaf already exist, directly return it.
     */
    const LeafNode* Find(int x, int y, int z, int depth_mask,
                         BranchNode* branch) const {
        int child = GetChildIndex(x, y, z, depth_mask);
        Node* node = branch->child_nodes_[child];
        if (node == NULL) return NULL;

        if (node->get_type() == BRANCH_NODE) {
            return Find(x, y, z, depth_mask >> 1,
                        reinterpret_cast<BranchNode*>(node));
        } else {
            return reinterpret_cast<LeafNode*>(node);
        }
    }

    /**
     * Erase a leaf node, return false if the node does not exist.
     */
    bool Erase(int x, int y, int z, int depth_mask, BranchNode* branch) {
        int child = GetChildIndex(x, y, z, depth_mask);
        Node* node = branch->child_nodes_[child];

        if (node == NULL) {
            return false;
        } else {
            switch (node->get_type()) {
            case BRANCH_NODE:
            {
                BranchNode* b = reinterpret_cast<BranchNode*>(node);
                if (Erase(x, y, z, depth_mask >> 1, b)) {
                    for (int i = 0; i < 8; ++i) {
                        if (b->child_nodes_[i])
                            return true;
                    }
                    branch_pool_.Free(b);
                    branch->child_nodes_[child] = NULL;
                    return true;
                } else {
                    return false;
                }
            }
            case LEAF_NODE:
            {
                --size_;
                LeafNode* leaf = reinterpret_cast<LeafNode*>(node);
                leaf_pool_.Free(leaf);
                branch->child_nodes_[child] = NULL;
                return true;
            }
            }
        }

        // Code never go here.
        assert(false && "Unreachable code.");
        return false;
    }

    /**
     * Initialize the octree.
     */
    void Initialize() {
        size_ = 0;
        root_ = branch_pool_.Allocate();

        SetDepth();
    }

    /**
     * Set the depth of octree.
     */
    void SetDepth() {
        // Find maximum amount of keys.
        int max_voxels = std::max(size1_, std::max(size2_, size3_));

        // Tree depth = amount of bits of max_voexls.
        depth_ = max_voxels == 0 ? 0
                                 : bits::Log2Floor(max_voxels) + 1;
    }

    /**
     * Get child node index using depth mask.
     */
    static int GetChildIndex(int x, int y, int z, int depth_mask) {
        return (((!!(x & depth_mask)) << 2) | ((!!(y & depth_mask)) << 1) |
                (!!(z & depth_mask)));
    }

    /**
     * @return the depth mask for fast searching the child.
     */
    int DepthMask() const {
        return 1 << (depth_ - 1);
    }

    int depth_;                          // Octree depth.
    BranchNode* root_;                   // The root node of octree.
    int size_;                           // The number of leaves.
    int size1_;                          // Size of the first dimension.
    int size2_;                          // Size of the second dimension.
    int size3_;                          // Size of the third dimension.
    ObjectPool<BranchNode> branch_pool_; // Memory pool for branch nodes.
    ObjectPool<LeafNode> leaf_pool_;     // Memory pool for leaf nodes.

private:
    DISALLOW_COPY_AND_ASSIGN(Octree);
};

} // namespace cl

#endif // UTIL_TREE_OCTREE_H_
