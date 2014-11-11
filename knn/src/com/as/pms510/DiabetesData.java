package com.as.pms510;

import java.util.Arrays;
import java.util.List;

public class DiabetesData {
	private final List<String> data;

	DiabetesData(List<String> d) {
		data = d;
	}

	public List<String> getData() {
		return data;
	}

	public static DiabetesData makeData(String csv_line) {
		final String[] d = csv_line.split(",");
		return new DiabetesData(Arrays.asList(d));
	}
}
