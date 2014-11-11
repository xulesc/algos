/*
The MIT License (MIT)
[OSI Approved License]
The MIT License (MIT)
Copyright (c) 2014 Daniel Glasson
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Modified from https://github.com/AReallyGoodName/OfflineReverseGeocode/tree/master/src/geocode/kdtree
*/

package com.as.pms510;

import java.util.Collections;
import java.util.List;

public class KDTree<T extends KComparator<T>> {
	private final KDNode<T> root;
	private int k;

	public KDTree(List<T> items, int k) {
		root = create(items, 0);
		this.k = k;
	}

	@SuppressWarnings("unchecked")
	public Object[] findKNearest(T search) {
		final Object[] nearest = new Object[k];
		final Object[] ret = new Object[k];
		// System.out.println(nearest.length);
		for (int i = 0; i < k; i++)
			nearest[i] = null;
		// System.out.println("processing: " + search);
		for (int i = 0; i < k; i++) {
			// System.out.println(">>>");
			nearest[i] = findNearest(root, search, 0, nearest);
			final T node = ((KDNode<T>) nearest[i]).nodeLocation;
			ret[i] = node;
			// System.out.println(">>>" + k + ":" + i + ":" + node);
		}
		return ret;
	}

	public T findNearest(T search) {
		return findNearest(root, search, 0, null).nodeLocation;
	}

	private KDNode<T> create(List<T> items, int depth) {
		if (items.isEmpty()) {
			return null;
		}
		Collections.sort(items, items.get(0).getComparator(depth % 3));
		final int currentIndex = items.size() / 2;
		return new KDNode<T>(
				create(items.subList(0, currentIndex), depth + 1),
				create(items.subList(currentIndex + 1, items.size()), depth + 1),
				items.get(currentIndex));
	}

	@SuppressWarnings("unchecked")
	private boolean check(Object[] nodes, KDNode<T> node) {
		if (nodes == null)
			return false;
		for (int i = 0; i < nodes.length; i++) {
			if (nodes[i] != null && node.equals((KDNode<T>) nodes[i])) {
				// System.out.println("found match: " + node + "==" +
				// (KDNode<T>)nodes[i]);
				return true;
			}
		}
		return false;
	}

	private KDNode<T> findNearest(KDNode<T> currentNode, T search, int depth,
			Object[] nearest) {
		final int direction = search.getComparator(depth % 3).compare(search,
				currentNode.nodeLocation);
		final KDNode<T> next = (direction < 0) ? currentNode.leftTree
				: currentNode.rightTree;
		final KDNode<T> other = (direction < 0) ? currentNode.rightTree
				: currentNode.leftTree;
		KDNode<T> best = (next == null) ? currentNode : findNearest(next,
				search, depth + 1, nearest);
		if (!check(nearest, currentNode)
				&& currentNode.nodeLocation.sqDistance(search) < best.nodeLocation
						.sqDistance(search)) {
			best = currentNode;
		}
		if (other == null)
			return best;

		if (!check(nearest, currentNode)
				&& currentNode.nodeLocation.axisSqDistance(search, depth % 3) < best.nodeLocation
						.sqDistance(search)) {
			final KDNode<T> possibleBest = findNearest(other, search,
					depth + 1, nearest);
			if (!check(nearest, possibleBest)
					&& possibleBest.nodeLocation.sqDistance(search) < best.nodeLocation
							.sqDistance(search)) {
				best = possibleBest;
			}
		}
		return best;
	}
}
