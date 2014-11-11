package com.as.pms510;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

class SimpleKNN implements NN {
	private final int k;
	private float[][] dX;
	private float[] dY;

	public SimpleKNN(int k) {
		this.k = k;
	}

	@Override
	public void fit(float[][] X, float[] Y) {
		dX = X;
		dY = Y;
	}

	@Override
	public float[] predict(float[][] X, final ExecutorService executor) {
		final float[] Y = new float[X.length];
		for (int i = 0; i < X.length; i++) {
			Y[i] = predict(X[i], executor);
		}
		return Y;
	}

	private float predict(final float[] x, final ExecutorService executor) {
		// make distances
		@SuppressWarnings("unused")
		final long t1 = System.currentTimeMillis();
		final double[] distances = new double[dX.length];
		// distribute dist calculation filling individual lists
		final Set<Callable<Pair>> callables = new HashSet<Callable<Pair>>();
		for (int i = 0; i < dX.length; i++) {
			final int index = i;
			callables.add(new Callable<Pair>() {
				@Override
				public Pair call() throws Exception {
					return new Pair(index, Utils.euclid_dist(x, dX[index]));
				}
			});
		}
		try {
			final List<Future<Pair>> futures = executor.invokeAll(callables);
			for (Future<Pair> future : futures) {
				distances[future.get().first] = future.get().second;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		// System.out.println("Distance calculation took "
		// + (System.currentTimeMillis() - t1) / 1000. + " seconds");
		// find top k classes
		@SuppressWarnings("unused")
		final long t2 = System.currentTimeMillis();
		final List<Float> klasses = new ArrayList<Float>();
		for (int i = 0; i < k; i++) {
			int index = -1;
			double min_dist = Double.MAX_VALUE;
			for (int j = 0; j < distances.length; j++) {
				if (min_dist > distances[j]) {
					index = j;
					min_dist = distances[j];
				}
			}
			klasses.add(dY[index]);
			distances[index] = Double.MAX_VALUE;
		}
		return Utils.getTopKlass(klasses);
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