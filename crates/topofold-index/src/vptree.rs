//! Vantage-Point Tree (VP-Tree) metric indexing.
//!
//! A VP-Tree indexes objects in an arbitrary metric space $(X, d)$ using spherical
//! decomposition. Search queries achieve sub-linear $O(\log N)$ expected time complexity
//! by using the triangle inequality to prune entire metric sub-trees.
//!
//! Reference:
//! - Yianilos, P. N. (1993). Data structures and algorithms for nearest neighbor search
//!   in general metric spaces. *SODA '93*, 311–321.

#![forbid(unsafe_code)]

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::metric::Metric;

/// A node in the Vantage-Point Tree.
#[derive(Debug, Clone)]
pub struct VpNode<T> {
    /// Vantage point item.
    pub item: T,
    /// Threshold radius (median distance to other elements in the sub-tree).
    pub threshold: f64,
    /// Sub-tree containing elements with distance $\le$ threshold.
    pub inside: Option<Box<VpNode<T>>>,
    /// Sub-tree containing elements with distance $>$ threshold.
    pub outside: Option<Box<VpNode<T>>>,
}

/// A Vantage-Point Tree index.
#[derive(Debug, Clone)]
pub struct VpTree<T, M> {
    root: Option<Box<VpNode<T>>>,
    metric: M,
    len: usize,
}

impl<T, M> VpTree<T, M>
where
    M: Metric<T>,
{
    /// Builds a new `VpTree` from a vector of items using the provided metric.
    #[must_use]
    pub fn new(items: Vec<T>, metric: M) -> Self {
        let len = items.len();
        let root = Self::build_recursive(items, &metric);
        Self { root, metric, len }
    }

    /// Number of items stored in the tree.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the tree contains no items.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Reference to the underlying metric.
    #[inline]
    #[must_use]
    pub fn metric(&self) -> &M {
        &self.metric
    }

    fn build_recursive(mut items: Vec<T>, metric: &M) -> Option<Box<VpNode<T>>> {
        if items.is_empty() {
            return None;
        }

        if items.len() == 1 {
            let item = items.pop().unwrap();
            return Some(Box::new(VpNode {
                item,
                threshold: 0.0,
                inside: None,
                outside: None,
            }));
        }

        let vp_item = items.swap_remove(0);

        // Partition remaining items around median distance to vp_item
        let median_idx = items.len() / 2;
        items.select_nth_unstable_by(median_idx, |a, b| {
            let d_a = metric.distance(&vp_item, a);
            let d_b = metric.distance(&vp_item, b);
            d_a.partial_cmp(&d_b).unwrap_or(Ordering::Equal)
        });

        let threshold = metric.distance(&vp_item, &items[median_idx]);

        // Split in-place into inside and outside partitions
        let outside_items = items.split_off(median_idx);
        let inside_items = items;

        let inside = Self::build_recursive(inside_items, metric);
        let outside = Self::build_recursive(outside_items, metric);

        Some(Box::new(VpNode {
            item: vp_item,
            threshold,
            inside,
            outside,
        }))
    }

    /// Searches for all items within `radius` distance of `query`.
    ///
    /// Prunes sub-trees using the metric triangle inequality:
    /// - Inside child explored if $d(q, v) - \text{radius} \le \text{threshold}$.
    /// - Outside child explored if $d(q, v) + \text{radius} > \text{threshold}$.
    #[must_use]
    pub fn range_search(&self, query: &T, radius: f64) -> Vec<(&T, f64)> {
        let mut results = Vec::new();
        if let Some(ref root) = self.root {
            Self::search_range_node(root, query, radius, &self.metric, &mut results);
        }
        results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
        results
    }

    fn search_range_node<'a>(
        node: &'a VpNode<T>,
        query: &T,
        radius: f64,
        metric: &M,
        results: &mut Vec<(&'a T, f64)>,
    ) {
        let dist = metric.distance(&node.item, query);

        if dist <= radius {
            results.push((&node.item, dist));
        }

        // Triangle inequality pruning:
        // Inside sphere: d(v, x) <= threshold
        // If dist - radius <= threshold, we must check inside
        if dist - radius <= node.threshold {
            if let Some(ref inside) = node.inside {
                Self::search_range_node(inside, query, radius, metric, results);
            }
        }

        // Outside shell: d(v, x) > threshold
        // If dist + radius > threshold, we must check outside
        if dist + radius > node.threshold {
            if let Some(ref outside) = node.outside {
                Self::search_range_node(outside, query, radius, metric, results);
            }
        }
    }

    /// Finds the $k$-nearest neighbors to `query`.
    #[must_use]
    pub fn k_nearest(&self, query: &T, k: usize) -> Vec<(&T, f64)> {
        if k == 0 || self.is_empty() {
            return Vec::new();
        }

        let mut heap: BinaryHeap<Neighbor<T>> = BinaryHeap::with_capacity(k);
        let mut tau = f64::INFINITY;

        if let Some(ref root) = self.root {
            Self::search_knn_node(root, query, k, &self.metric, &mut heap, &mut tau);
        }

        let mut results = Vec::with_capacity(heap.len());
        while let Some(n) = heap.pop() {
            results.push((n.item, n.dist));
        }
        results.reverse(); // Smallest distance first
        results
    }

    fn search_knn_node<'a>(
        node: &'a VpNode<T>,
        query: &T,
        k: usize,
        metric: &M,
        heap: &mut BinaryHeap<Neighbor<'a, T>>,
        tau: &mut f64,
    ) {
        let dist = metric.distance(&node.item, query);

        if dist < *tau {
            if heap.len() == k {
                heap.pop();
            }
            heap.push(Neighbor {
                dist,
                item: &node.item,
            });

            if heap.len() == k {
                *tau = heap.peek().map_or(f64::INFINITY, |top| top.dist);
            }
        }

        // Determine which branch to search first based on query position
        let search_inside_first = dist < node.threshold;

        if search_inside_first {
            if dist - *tau <= node.threshold {
                if let Some(ref inside) = node.inside {
                    Self::search_knn_node(inside, query, k, metric, heap, tau);
                }
            }
            if dist + *tau > node.threshold {
                if let Some(ref outside) = node.outside {
                    Self::search_knn_node(outside, query, k, metric, heap, tau);
                }
            }
        } else {
            if dist + *tau > node.threshold {
                if let Some(ref outside) = node.outside {
                    Self::search_knn_node(outside, query, k, metric, heap, tau);
                }
            }
            if dist - *tau <= node.threshold {
                if let Some(ref inside) = node.inside {
                    Self::search_knn_node(inside, query, k, metric, heap, tau);
                }
            }
        }
    }
}

/// Helper struct for $k$-NN max-heap ordered by distance.
struct Neighbor<'a, T> {
    dist: f64,
    item: &'a T,
}

impl<'a, T> PartialEq for Neighbor<'a, T> {
    fn eq(&self, other: &Self) -> bool {
        self.dist == other.dist
    }
}

impl<'a, T> Eq for Neighbor<'a, T> {}

impl<'a, T> PartialOrd for Neighbor<'a, T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a, T> Ord for Neighbor<'a, T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.dist
            .partial_cmp(&other.dist)
            .unwrap_or(Ordering::Equal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Euclidean1D;
    impl Metric<f64> for Euclidean1D {
        fn distance(&self, a: &f64, b: &f64) -> f64 {
            (a - b).abs()
        }
    }

    #[test]
    fn test_vptree_1d_search() {
        let points = vec![1.0, 3.0, 7.0, 8.0, 12.0, 15.0, 20.0];
        let tree = VpTree::new(points, Euclidean1D);

        // Range query at 7.5 with radius 1.0: should return 7.0 and 8.0
        let range_hits = tree.range_search(&7.5, 1.0);
        let hit_vals: Vec<f64> = range_hits.into_iter().map(|(&x, _)| x).collect();
        assert_eq!(hit_vals, vec![7.0, 8.0]);

        // 3-NN query around 10.0: should return 8.0 (dist 2.0), 12.0 (dist 2.0), 7.0 (dist 3.0)
        let knn_hits = tree.k_nearest(&10.0, 3);
        assert_eq!(knn_hits.len(), 3);
        let knn_vals: Vec<f64> = knn_hits.into_iter().map(|(&x, _)| x).collect();
        assert!(knn_vals.contains(&8.0));
        assert!(knn_vals.contains(&12.0));
        assert!(knn_vals.contains(&7.0));
    }
}
