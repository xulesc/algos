package com.as.pms510;

import java.util.ArrayList;
import java.util.List;

public class DiabetesData2 {
	public static int NUMERICS = 0;
	public static int NOMINALS = 1;
	public static int ALL = 2;
	private List<Float> numerics;
	private List<Float> nominals;
	private float klass;

	DiabetesData2(List<String> s_numerics, List<String> s_nominals,
			String s_klass) {
		numerics = new ArrayList<Float>();
		nominals = new ArrayList<Float>();
		//
		for (String numeric : s_numerics)
			numerics.add(Float.parseFloat(numeric));
		for (String nominal : s_nominals)
			nominals.add(Float.parseFloat(nominal));
		klass = Float.parseFloat(s_klass);
	}

	public static float[][] X(List<DiabetesData2> data, int type) {
		int dim_count = 0;
		if (type == NUMERICS || type == ALL)
			dim_count += data.get(0).numerics.size();
		if (type == NOMINALS || type == ALL)
			dim_count += data.get(0).nominals.size();

		final float[][] x = new float[data.size()][dim_count];
		int row = 0;
		for (DiabetesData2 d : data) {
			int col = 0;
			if (type == NUMERICS || type == ALL) {
				for (float f : d.numerics)
					x[row][col++] = f;
			}
			if (type == NOMINALS || type == ALL) {
				for (float f : d.nominals)
					x[row][col++] = f;
			}
			row += 1;
		}
		return x;
	}

	public static float[] Y(List<DiabetesData2> data) {
		final float[] ret = new float[data.size()];
		for (int i = 0; i < data.size(); i++) {
			final DiabetesData2 abaloneData = data.get(i);
			ret[i] = abaloneData.klass;
		}
		return ret;
	}

	public static DiabetesData2 makeData(List<String> s_numerics,
			List<String> s_nominals, String s_klass) {
		for (String numeric : s_numerics)
			if (numeric.equals("?") || numeric.startsWith("V")
					|| numeric.startsWith("E"))
				return null;
		//
		return new DiabetesData2(s_numerics, s_nominals, s_klass);
	}
}
