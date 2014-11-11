package com.as.pms510;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class Utils {
	public static double euclid_dist(float[] x1, float[] x2) {
		float dist = 0;
		for (int i = 0; i < x1.length; i++) {
			dist += Math.pow((x1[i] - x2[i]), 2);
		}
		return Math.sqrt(dist);
	}

	public static int find_mean(float[] x, float[][] means) {
		final double[] distances = new double[means.length];
		for (int i = 0; i < means.length; i++)
			distances[i] = Utils.euclid_dist(x, means[i]);
		int index = -1;
		double min_dist = Double.MAX_VALUE;
		for (int i = 0; i < distances.length; i++) {
			if (distances[i] < min_dist) {
				index = i;
				min_dist = distances[i];
			}
		}
		return index;
	}

	public static float getTopKlass(List<Float> klasses) {
		final Map<Float, Integer> counts = new HashMap<Float, Integer>();
		for (int i = 0; i < klasses.size(); i++) {
			Integer v = counts.get(klasses.get(i));
			if (v == null) {
				v = 0;
			}
			counts.put(klasses.get(i), v + 1);
		}
		float klass = Float.MIN_VALUE;
		int max_count = 0;
		for (Entry<Float, Integer> entry : counts.entrySet()) {
			if (entry.getValue() > max_count) {
				klass = entry.getKey();
				max_count = entry.getValue();
			}
		}
		return klass;
	}

}
