package com.as.pms510;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class DiabetesML {
	/**
	 * apply transformation described in paper:
	 * 
	 * 'Beata Strack, Jonathan P. DeShazo, Chris Gennings, Juan L. Olmo,
	 * Sebastian Ventura, Krzysztof J. Cios, and John N. Clore, “Impact of HbA1c
	 * Measurement on Hospital Readmission Rates: Analysis of 70,000 Clinical
	 * Database Patient Records,” BioMed Research International, vol. 2014,
	 * Article ID 781670, 11 pages, 2014.'
	 * 
	 * in addition separate the features in numeric and string types.
	 */
	private final static String runModePreProcess = "prep";
	/**
	 * convert columns to nominal values for the file containing data dimensions
	 * that is expected to be strings or nominal values
	 */
	private final static String runModeNominalize = "nomn";
	private final static String runModeClusterize = "clust";
	private final static String runModeClusterize2 = "clust2";
	private final static String runModeBenchmarks = "bench";
	private final static String runModeReductions = "reduce";
	/**
	 * column names to be excluded from inclusion into the pre-processed data
	 * file generated as per the above paper due to lack on information.
	 */
	private final static Set<String> excludeColumnNames = new HashSet<String>(
			Arrays.asList("weight", "payer_code", "medical_specialty"));
	/**
	 * discharge codes that refer to the patient being discharged to Hospice
	 */
	private final static Set<String> excludedDischarge = new HashSet<String>(
			Arrays.asList("11", "13", "14", "19", "20", "21"));
	/**
	 * columns that are already nominalized in the original dataset
	 */
	private final static Set<String> forceNominal = new HashSet<String>(
			Arrays.asList("discharge_disposition_id", "admission_source_id",
					"admission_type_id"));

	private final static String MAPPINGS_EXT = ".mappings";
	private final static String NUMERICS_EXT = ".num";
	private final static String NON_NUMERICS_EXT = ".non_num";
	private final static String NOMINALS_EXT = NON_NUMERICS_EXT + ".nominal";
	private final static String KLASS_EXT = ".klass";

	public static void printUsage() {
		System.out.println("usage: ");
	}

	public static void main(String args[]) throws IOException {
		if (args.length < 2) {
			printUsage();
			System.exit(-1);
		}

		final String dataFileName = args[0];
		final String runMode = args[1];

		if (runMode.equals(runModePreProcess)) {
			new DiabetesML().preProcess(dataFileName);
		} else if (runMode.equals(runModeNominalize)) {
			new DiabetesML().nominalize(dataFileName);
		} else if (runMode.equals(runModeClusterize)) {
			new DiabetesML().clusterize(dataFileName);
		} else if (runMode.equals(runModeBenchmarks)) {
			new DiabetesML().benchmark(dataFileName);
		} else if (runMode.equals(runModeReductions)) {
			new DiabetesML().reduce(dataFileName);
		} else if (runMode.equals(runModeClusterize2)) {
			new DiabetesML().clusterize2(dataFileName);
		} else {
			System.out.println("nothing to do");
			System.exit(-2);
		}
	}

	private List<String> readKlasses(String dataFileName) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(dataFileName
				+ KLASS_EXT));
		String strLine = "";
		final List<String> d_klasses = new ArrayList<String>();
		int row = 0;
		while ((strLine = br.readLine()) != null) {
			if (row++ == 0)
				continue;
			d_klasses.add(strLine);
		}
		br.close();
		return d_klasses;
	}

	private List<List<String>> readNumerics(String dataFileName)
			throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(dataFileName
				+ NUMERICS_EXT));
		String strLine = "";
		int row = 0;
		final List<List<String>> d_numerics = new ArrayList<List<String>>();
		while ((strLine = br.readLine()) != null) {
			if (row++ == 0)
				continue;
			d_numerics.add(Arrays.asList(strLine.split(",")));
		}
		br.close();
		return d_numerics;
	}

	private List<List<String>> readNominals(String dataFileName)
			throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(dataFileName
				+ NOMINALS_EXT));
		String strLine = "";
		int row = 0;
		final List<List<String>> d_nominals = new ArrayList<List<String>>();
		while ((strLine = br.readLine()) != null) {
			if (row++ == 0)
				continue;
			d_nominals.add(Arrays.asList(strLine.split(",")));
		}
		br.close();
		return d_nominals;
	}

	private void benchmark(String dataFileName) throws IOException {
		final List<DiabetesData2> data = loadData(dataFileName);
		final List<DiabetesData2> benchmark_data = data.subList(0,
				data.size() / 10);
		final Pair trg_tst_sets = getTrainingTestSet(benchmark_data);
		System.out.println("#trg set size: " + trg_tst_sets.first.size());
		System.out.println("#trg set size: " + trg_tst_sets.second.size());
		//
		final int core_count = Runtime.getRuntime().availableProcessors();
		benchmark_cc(core_count, trg_tst_sets.first, trg_tst_sets.second);
		benchmark_k(core_count, trg_tst_sets.first, trg_tst_sets.second);
		benchmark_fset(core_count, trg_tst_sets.first, trg_tst_sets.second);
	}

	private Pair2 loadClust2TrgData(String dataFileName) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(dataFileName
				+ ".reduced"));
		String strLine = "";
		final List<List<Float>> data = new ArrayList<List<Float>>();
		while ((strLine = br.readLine()) != null) {
			final List<Float> row = new ArrayList<Float>();
			for (String s : strLine.split(","))
				row.add(Float.parseFloat(s));
			data.add(row);
		}
		br.close();
		final int dimC = data.get(0).size() - 1;
		final float[][] x = new float[data.size()][dimC];
		final float[] y = new float[data.size()];
		for (int i = 0; i < data.size(); i++) {
			y[i] = data.get(i).get(dimC);
			for (int j = 0; j < dimC; j++)
				x[i][j] = data.get(i).get(j);
		}
		return new Pair2(x, y);
	}

	private void clusterize2(String dataFileName) throws IOException {
		//
		final Pair2 trgPair = loadClust2TrgData(dataFileName);
		//
		final List<DiabetesData2> data = loadData(dataFileName);
		final Pair pair = getTrainingTestSet(data);
		final List<DiabetesData2> tst = pair.second;
		System.out.println("#tsts set size: " + tst.size());
		float[][] tst_x = DiabetesData2.X(tst, DiabetesData2.ALL);
		float[] tst_y = DiabetesData2.Y(tst);
		System.out.println("#tst dim count: " + tst_x[0].length);
		//
		normalize(tst_x);
		//
		final int core_count = Runtime.getRuntime().availableProcessors();
		final ExecutorService executor = Executors
				.newFixedThreadPool(core_count);
		final long t1 = System.currentTimeMillis();
		final int k = 1;
		final NN knn = new KDTreeNN(k);
		final float[][] trg_x = trgPair.first;
		System.out.println("#trg dim count: " + trg_x[0].length);
		knn.fit(trg_x, trgPair.second);
		final float[] predicted = knn.predict(tst_x, executor);
		score(predicted, tst_y, DiabetesData2.ALL, core_count, true, k, t1,
				"clust2");
		executor.shutdown();
		//
		float[][] x_tst_x = tst_x; // DiabetesData2.X(tst, DiabetesData2.ALL);
		final PrintWriter writer = new PrintWriter(dataFileName
				+ ".reduced.clust2", "UTF-8");
		for (int i = 0; i < x_tst_x.length; i++) {
			final List<Float> o = new ArrayList<Float>();
			o.add(predicted[i]);
			o.add(tst_y[i]);
			for (int j = 0; j < x_tst_x[i].length; j++)
				o.add(x_tst_x[i][j]);
			writer.println(join(o, " "));
		}
		writer.close();
		//
	}

	private void clusterize(String dataFileName) throws IOException {
		final List<DiabetesData2> data = loadData(dataFileName);
		final Pair pair = getTrainingTestSet(data);
		final List<DiabetesData2> trg = pair.first;
		final List<DiabetesData2> tst = pair.second;
		System.out.println("#trg set size: " + trg.size());
		System.out.println("#trg set size: " + tst.size());
		//
		float[][] trg_x = DiabetesData2.X(trg, DiabetesData2.NUMERICS);
		float[] trg_y = DiabetesData2.Y(trg);
		float[][] tst_x = DiabetesData2.X(tst, DiabetesData2.NUMERICS);
		float[] tst_y = DiabetesData2.Y(tst);
		System.out.println("#dim count: " + trg_x[0].length);
		//
		normalize(trg_x);
		normalize(tst_x);
		//
		final int core_count = Runtime.getRuntime().availableProcessors();
		final ExecutorService executor = Executors
				.newFixedThreadPool(core_count);
		final long t1 = System.currentTimeMillis();
		final NN knn = new KDTreeNN(1);
		knn.fit(trg_x, trg_y);
		final float[] predicted = knn.predict(tst_x, executor);
		score(predicted, tst_y, DiabetesData2.NUMERICS, core_count, true, 1,
				t1, "clust");
		executor.shutdown();
		//
		float[][] x_tst_x = DiabetesData2.X(tst, DiabetesData2.NOMINALS);
		final PrintWriter writer = new PrintWriter(dataFileName + ".clust2",
				"UTF-8");
		for (int i = 0; i < x_tst_x.length; i++) {
			final List<Float> o = new ArrayList<Float>();
			o.add(predicted[i]);
			o.add(tst_y[i]);
			for (int j = 0; j < x_tst_x[i].length; j++)
				o.add(x_tst_x[i][j]);
			writer.println(join(o, " "));
		}
		writer.close();
	}

	private void reduce(String dataFileName) throws IOException {
		final List<DiabetesData2> data = loadData(dataFileName);
		float[][] trg_x = DiabetesData2.X(data, DiabetesData2.ALL);
		float[] trg_y = DiabetesData2.Y(data);
		normalize(trg_x);
		final KDTreeNN knn = new KDTreeNN(1);
		final List<float[]> new_data = new ArrayList<float[]>();
		final List<Float> new_data_klass = new ArrayList<Float>();
		final PrintWriter writer = new PrintWriter(dataFileName + ".reduced",
				"UTF-8");
		for (int i = 0; i < trg_x.length; i++) {
			if (i % 1000 == 0)
				System.out.println("at index: " + i);
			if (i != 0 && knn.nearestPoint(trg_x[i]) == trg_y[i])
				continue;
			new_data.add(trg_x[i]);
			new_data_klass.add(trg_y[i]);
			//
			final float[][] tr_x = new float[new_data.size()][trg_x[0].length];
			final float[] tr_y = new float[new_data.size()];
			for (int j = 0; j < new_data.size(); j++)
				for (int k = 0; k < trg_x[0].length; k++)
					tr_x[j][k] = new_data.get(j)[k];
			for (int j = 0; j < new_data.size(); j++)
				tr_y[j] = new_data_klass.get(j);
			knn.fit(tr_x, tr_y);
			//
			final List<Float> o = new ArrayList<Float>();
			for (float f : trg_x[i])
				o.add(f);
			o.add(trg_y[i]);
			writer.println(join(o, ","));
		}
		//
		writer.close();
	}

	private List<DiabetesData2> loadData(String dataFileName)
			throws IOException {
		// read klasses, numerics and nominals
		final List<String> d_klasses = readKlasses(dataFileName);
		final List<List<String>> d_numerics = readNumerics(dataFileName);
		final List<List<String>> d_nominals = readNominals(dataFileName);
		//
		final List<DiabetesData2> data = new ArrayList<DiabetesData2>();
		for (int i = 0; i < d_klasses.size(); i++) {
			final DiabetesData2 data_point = DiabetesData2.makeData(
					d_numerics.get(i), d_nominals.get(i), d_klasses.get(i));
			if (data_point == null)
				continue;
			data.add(data_point);
		}
		System.out.println("#dataset size: " + data.size());
		//
		Collections.shuffle(data);
		return data;
	}

	private void benchmark_cc(int core_count, List<DiabetesData2> trg,
			List<DiabetesData2> tst) {
		// test CC for kdtree
		for (int cc = 1; cc <= 2 * core_count; cc++) {
			final ExecutorService executor = Executors.newFixedThreadPool(cc);
			classify(trg, tst, DiabetesData2.NUMERICS, executor, true, cc, 1,
					"kdtree-cc");
			executor.shutdown();
		}
		// test CC for non-kdtree
		for (int cc = 1; cc <= 2 * core_count; cc++) {
			final ExecutorService executor = Executors.newFixedThreadPool(cc);
			classify(trg, tst, DiabetesData2.NUMERICS, executor, false, cc, 1,
					"non-kdtree-cc");
			executor.shutdown();
		}
	}

	private void benchmark_k(int core_count, List<DiabetesData2> trg,
			List<DiabetesData2> tst) {
		final ExecutorService executor = Executors
				.newFixedThreadPool(core_count);
		for (int i = 2; i <= 15; i++) {
			classify(trg, tst, DiabetesData2.NUMERICS, executor, true,
					core_count, i, "kdtree-k");
			// classify(trg, tst, DiabetesData2.NOMINALS, executor, true,
			// core_count, i);
			// classify(trg, tst, DiabetesData2.ALL, executor, true, core_count,
			// i);
		}
		executor.shutdown();
	}

	private void benchmark_fset(int core_count, List<DiabetesData2> trg,
			List<DiabetesData2> tst) {
		final ExecutorService executor = Executors
				.newFixedThreadPool(core_count);
		classify(trg, tst, DiabetesData2.NUMERICS, executor, true, core_count,
				1, "kdtree-numerics");
		classify(trg, tst, DiabetesData2.NOMINALS, executor, true, core_count,
				1, "kdtree-nominals");
		classify(trg, tst, DiabetesData2.ALL, executor, true, core_count, 1,
				"kdtree-all");
		executor.shutdown();
	}

	private Pair getTrainingTestSet(List<DiabetesData2> data) {
		final int trg_size = (int) Math.floor(data.size() * 95 / 100.0);
		final List<DiabetesData2> trg = data.subList(0, trg_size);
		final List<DiabetesData2> tst = data.subList(trg_size, data.size());
		return new Pair(trg, tst);
	}

	private float[][] normalize(float[][] d) {
		for (int col = 0; col < d[0].length; col++) {
			float max = -1;
			for (int row = 0; row < d.length; row++) {
				if (d[row][col] > max)
					max = d[row][col];
			}
			if (max != 0) {
				for (int row = 0; row < d.length; row++) {
					d[row][col] /= max;
				}
			}
		}
		return d;
	}

	private void classify(final List<DiabetesData2> trg_data,
			final List<DiabetesData2> tst_data, int type,
			final ExecutorService executor, final boolean use_kdtree, int cc,
			int k, String logId) {
		float[][] trg_x = DiabetesData2.X(trg_data, type);
		float[] trg_y = DiabetesData2.Y(trg_data);
		float[][] tst_x = DiabetesData2.X(tst_data, type);
		float[] tst_y = DiabetesData2.Y(tst_data);
		//
		normalize(trg_x);
		normalize(tst_x);
		//
		classify(k, trg_x, tst_x, trg_y, tst_y, type, executor, use_kdtree, cc,
				logId);
	}

	private void classify(int k, float[][] trg_x, float[][] tst_x,
			float[] trg_y, float[] tst_y, int type,
			final ExecutorService executor, boolean use_kdtree, int cc,
			String logId) {
		final long t1 = System.currentTimeMillis();
		final NN knn = (use_kdtree) ? new KDTreeNN(k) : new SimpleKNN(k);
		knn.fit(trg_x, trg_y);
		final float[] predicted = knn.predict(tst_x, executor);
		score(predicted, tst_y, type, cc, use_kdtree, k, t1, logId);
	}

	private static void score(float[] predicted_y, float[] actual_y, int type,
			int cc, boolean use_kdtree, int k, long t1, String logId) {
		int correct = 0;
		for (int i = 0; i < predicted_y.length; i++)
			if (predicted_y[i] == actual_y[i])
				correct++;
		int a = 0, b = 0, c = 0;
		for (int i = 0; i < predicted_y.length; i++) {
			if (predicted_y[i] == 1 && actual_y[i] == 1)
				a += 1;
			if (predicted_y[i] == 1 && actual_y[i] == 0)
				b += 1;
			if (predicted_y[i] == 0 && actual_y[i] == 1)
				c += 1;
		}
		System.out.println(a + ":" + b + ":" + c);
		final long tdiff = System.currentTimeMillis() - t1;
		System.out
				.println("#" + logId + " " + type + ", cc: " + cc
						+ ", kdtree: " + use_kdtree + ", k: " + k + ", score: "
						+ 1.0 * correct / predicted_y.length + ", Prec: " + 1.0
						* a / (a + b) + ", Rec: " + 1.0 * a / (a + c)
						+ ", TDiff: " + (tdiff / 1000.));
	}

	private void nominalize(String dataFileName) throws IOException {
		// read file and build unique value map
		BufferedReader br = new BufferedReader(new FileReader(dataFileName));
		String strLine = "";
		final List<Map<String, Integer>> value_mappings = new ArrayList<Map<String, Integer>>();
		while ((strLine = br.readLine()) != null) {
			final String[] cols_data = strLine.split(",");
			// first line is header. find out column count and initialize list.
			// dont use last column it is the klass
			if (value_mappings.isEmpty()) {
				for (int i = 0; i < cols_data.length - 1; i++) {
					value_mappings.add(new HashMap<String, Integer>());
				}
				continue;
			}
			// put all unique values per column into the map
			for (int i = 0; i < cols_data.length - 1; i++) {
				final String col_data = cols_data[i];
				if (!value_mappings.get(i).containsKey(col_data))
					value_mappings.get(i).put(col_data,
							value_mappings.get(i).size());
			}
		}
		br.close();
		// read file and write out values replaced by nominal values
		br = new BufferedReader(new FileReader(dataFileName));
		BufferedWriter bw = new BufferedWriter(new FileWriter(dataFileName
				+ ".nominal"));
		strLine = "";
		int lcount = 0;
		while ((strLine = br.readLine()) != null) {
			if (++lcount == 1) {
				bw.write(strLine + "\n");
				continue;
			}
			final String[] data = strLine.split(",");
			final List<Integer> output = new ArrayList<Integer>();
			for (int i = 0; i < data.length - 1; i++)
				output.add(value_mappings.get(i).get(data[i]));
			bw.write(join(output, ",") + "\n");
		}
		br.close();
		bw.close();
		// write out the mappings
		bw = new BufferedWriter(new FileWriter(dataFileName + MAPPINGS_EXT));
		for (int i = 0; i < value_mappings.size(); i++) {
			bw.write(i + "------------------\n");
			for (Map.Entry<String, Integer> entry : value_mappings.get(i)
					.entrySet()) {
				bw.write("\t" + entry.getKey() + ":" + entry.getValue() + "\n");
			}
		}
		bw.close();
	}

	private static String join(List<?> list, String delim) {
		final StringBuilder sb = new StringBuilder();
		String loopDelim = "";
		for (Object s : list) {
			sb.append(loopDelim);
			sb.append(s);
			loopDelim = delim;
		}
		return sb.toString();
	}

	private void preProcess(String dataFileName) throws IOException {
		// read original diabetes dataset
		final BufferedReader br = new BufferedReader(new FileReader(
				dataFileName));
		String sCurrentLine;
		final Map<String, List<String>> data = new LinkedHashMap<String, List<String>>();
		int index = 0;
		final List<String> headers = new ArrayList<String>();
		while ((sCurrentLine = br.readLine()) != null) {
			if (index++ == 0) {
				headers.addAll(Arrays.asList(sCurrentLine.split(",")));
				for (String header : headers)
					data.put(header, new ArrayList<String>());
				continue;
			}
			final List<String> d = DiabetesData.makeData(sCurrentLine)
					.getData();
			for (int i = 0; i < d.size(); i++) {
				final List<String> value = data.get(headers.get(i));
				value.add(d.get(i));
			}
		}
		br.close();
		// exclude by columns
		for (final String exclude : excludeColumnNames)
			data.remove(exclude);
		// exclude by rows
		final List<List<String>> dataset = new ArrayList<List<String>>();
		final List<String> dataset_headers = new ArrayList<String>(
				data.keySet());
		final int row_count = data.get(headers.get(0)).size();
		final Set<Integer> seen_pids = new HashSet<Integer>();
		for (int i = 0; i < row_count; i++) {
			final List<String> row = new ArrayList<String>();
			for (String key : data.keySet()) {
				row.add(data.get(key).get(i));
			}
			final int pid = Integer.parseInt(row.get(dataset_headers
					.indexOf("patient_nbr")));
			if (seen_pids.contains(pid))
				continue;
			final String disch = row.get(dataset_headers
					.indexOf("discharge_disposition_id"));
			if (excludedDischarge.contains(disch))
				continue;
			seen_pids.add(pid);
			dataset.add(row);
		}
		System.out.println("#dataset size: " + dataset.size());
		System.out.println("#dataset dims: " + dataset_headers.size());
		// segregate
		final List<Integer> classification = new ArrayList<Integer>();
		final List<List<String>> numeric_dims = new ArrayList<List<String>>();
		final List<List<String>> non_numeric_dims = new ArrayList<List<String>>();
		for (List<String> row : dataset)
			classification.add((row.get(row.size() - 1).equals("<30")) ? 1 : 0);
		final List<Integer> numeric_columns = new ArrayList<Integer>();
		final List<Integer> exclude_columns = new ArrayList<Integer>();
		for (int i = 0; i < dataset.get(0).size(); i++) {
			if (dataset_headers.get(i).equals("readmitted")) {
				exclude_columns.add(i);
				continue;
			}
			if (dataset_headers.get(i).equals("encounter_id")) {
				exclude_columns.add(i);
				continue;
			}
			if (dataset_headers.get(i).equals("patient_nbr")) {
				exclude_columns.add(i);
				continue;
			}
			//
			if (forceNominal.contains(dataset_headers.get(i)))
				continue;
			try {
				Double.parseDouble(dataset.get(0).get(i));
			} catch (Exception e) {
				continue;
			}
			numeric_columns.add(i);
		}
		System.out.println("#numeric columns: " + numeric_columns.size());
		for (List<String> row : dataset) {
			final List<String> numeric_row = new ArrayList<String>();
			final List<String> non_numeric_row = new ArrayList<String>();
			for (int i = 0; i < row.size(); i++) {
				if (exclude_columns.contains(i))
					continue;
				if (numeric_columns.contains(i)) {
					numeric_row.add(row.get(i));
				} else {
					non_numeric_row.add(row.get(i));
				}
			}
			numeric_dims.add(numeric_row);
			non_numeric_dims.add(non_numeric_row);
		}
		// write out
		BufferedWriter bw = new BufferedWriter(new FileWriter(dataFileName
				+ KLASS_EXT));
		bw.write("klass\n");
		for (Integer klass : classification)
			bw.write(klass + "\n");
		bw.close();
		bw = new BufferedWriter(new FileWriter(dataFileName + NUMERICS_EXT));
		final List<String> num_headers = new ArrayList<String>();
		for (Integer numeric_column : numeric_columns)
			num_headers.add(dataset_headers.get(numeric_column));
		bw.write(join(num_headers, ",") + "\n");
		for (List<String> num_row : numeric_dims)
			bw.write(join(num_row, ",") + "\n");
		bw.close();
		bw = new BufferedWriter(new FileWriter(dataFileName + NON_NUMERICS_EXT));
		final List<String> non_num_headers = new ArrayList<String>();
		for (int i = 0; i < dataset_headers.size(); i++)
			if (!numeric_columns.contains(i) && !exclude_columns.contains(i))
				non_num_headers.add(dataset_headers.get(i));
		bw.write(join(non_num_headers, ",") + "\n");
		for (List<String> non_num_row : non_numeric_dims)
			bw.write(join(non_num_row, ",") + "\n");
		bw.close();
		//
	}

	class Pair {
		final public List<DiabetesData2> first;
		final public List<DiabetesData2> second;

		Pair(List<DiabetesData2> f, List<DiabetesData2> s) {
			first = f;
			second = s;
		}
	}

	class Pair2 {
		final public float[][] first;
		final public float[] second;

		Pair2(float[][] f, float[] s) {
			first = f;
			second = s;
		}
	}
}
