package com.as.pms510;

import java.util.concurrent.ExecutorService;

public interface NN {

	public abstract void fit(float[][] X, float[] Y);

	public abstract float[] predict(float[][] X, ExecutorService executor);

}