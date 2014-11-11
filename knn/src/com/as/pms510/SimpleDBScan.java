package com.as.pms510;

import java.util.ArrayList;
import java.util.List;

/*
 * Implementation based on: http://en.wikipedia.org/wiki/DBSCAN#Algorithm
 */
public class SimpleDBScan {
    private final float eps;
    private final int minPts;
    private List<List<Integer>> neighbourHoods;
    private boolean[] visited;
    private boolean[] clustered;
    private List<List<Integer>> clusters;
    private float[][] cluster_centers;

    public SimpleDBScan(float e, int m) {
        eps = e;
        minPts = m;
    }

    public void fit(float[][] X) {
        // TODO this can be a triangular matrix to save memory & time?
        neighbourHoods = new ArrayList<List<Integer>>();
        for (int i = 0; i < X.length; i++) {
            neighbourHoods.add(new ArrayList<Integer>());
            for (int j = 0; j < X.length; j++) {
                final double dist = Utils.euclid_dist(X[i], X[j]);
                if (dist < eps)
                    neighbourHoods.get(i).add(j);
            }
        }
        // init
        visited = new boolean[X.length];
        clustered = new boolean[X.length];
        for (int i = 0; i < X.length; i++) {
            visited[i] = false;
            clustered[i] = false;
        }
        clusters = new ArrayList<List<Integer>>();
        //
        doDBScan(X);
        // make cluster centers
        cluster_centers = new float[clusters.size()][X[0].length];
        int cindex = 0;
        for(List<Integer> cluster : clusters) {
            final float[] cluster_center = new float[X[0].length];
            for(int point : cluster) {
                final float[] x = X[point];
                for(int i = 0; i < x.length; i++)
                    cluster_center[i] += x[i];
            }
            for(int i = 0; i < cluster_center.length; i++) {
                cluster_centers[cindex][i] = cluster_center[i]/cluster.size();
            }
            cindex++;
        }
    }
    
    float[] predict(float[][] X) {
        final float[] predicted = new float[X.length];
        
        for(int i = 0; i < X.length; i++)
            predicted[i] = Utils.find_mean(X[i], cluster_centers);
        
        return predicted;                   
    }    

    private void doDBScan(float[][] X) {
        for (int i = 0; i < X.length; i++) {
            if (visited[i])
                continue;
            visited[i] = true;
            final List<Integer> neighbourPts = regionQuery(i);
            if (neighbourPts.size() < minPts) {
                continue;
            }
            final List<Integer> cluster = new ArrayList<Integer>();
            clusters.add(cluster);
            expandCluster(i, neighbourPts, cluster);
        }
    }

    private void expandCluster(int P, List<Integer> neighbourPts,
            List<Integer> C) {
        C.add(P);
        clustered[P] = true;
        int i = 0;
        while (i < neighbourPts.size()) {
            final int neighbourPt = neighbourPts.get(i++);
            if (!visited[neighbourPt]) {
                visited[neighbourPt] = true;
                final List<Integer> neighbourPts1 = regionQuery(neighbourPt);
                if (neighbourPts1.size() > minPts) {
                    neighbourPts.removeAll(neighbourPts1);
                    neighbourPts.addAll(neighbourPts1);
                }
            }
            if (!clustered[neighbourPt]) {
                clustered[neighbourPt] = true;
                C.add(neighbourPt);
            }
        }
    }

    private List<Integer> regionQuery(int P) {
        return new ArrayList<Integer>(neighbourHoods.get(P));
    }
}
