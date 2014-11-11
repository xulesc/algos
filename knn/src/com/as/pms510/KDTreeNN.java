package com.as.pms510;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

public class KDTreeNN implements NN {
	KDTree<DataKDNode> kdTree;
	final private int k;

	KDTreeNN(int k) {
		this.k = k;
	}

	@Override
	public void fit(float[][] X, float[] Y) {
		final List<DataKDNode> trg = new ArrayList<DataKDNode>();
		int index = 0;
		for (float[] point : X)
			trg.add(new DataKDNode(point, Y[index++]));
		kdTree = new KDTree<DataKDNode>(trg, k);
	}

	public float nearestPoint(float[] x) {
		final DataKDNode n = kdTree.findNearest(new DataKDNode(x, -1));
		if (n == null)
			return -5000;
		return n.getKlass();
	}

	@Override
	public float[] predict(float[][] X, final ExecutorService executor) {
		final List<DataKDNode> tst = new ArrayList<DataKDNode>();
		for (float[] point : X)
			tst.add(new DataKDNode(point, -1));

		final float[] predicted = new float[X.length];
		// distribute dist calculation filling individual lists

		final Set<Callable<Pair>> callables = new HashSet<Callable<Pair>>();
		for (int i = 0; i < tst.size(); i++) {
			final DataKDNode point = tst.get(i);
			final int index = i;
			callables.add(new Callable<Pair>() {
				@Override
				public Pair call() throws Exception {
					final Object[] neighs = kdTree.findKNearest(point);
					final List<Float> klasses = new ArrayList<Float>();
					for (Object n : neighs) {
						klasses.add(((DataKDNode) n).getKlass());
					}
					return new Pair(index, Utils.getTopKlass(klasses));
				}
			});
		}
		try {
			for (Future<Pair> future : executor.invokeAll(callables)) {
				predicted[future.get().first] = (float) future.get().second;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		//
		return predicted;
	}

	class Pair {
		public final int first;
		public final double second;

		Pair(int f, double s) {
			first = f;
			second = s;
		}
	}

}
