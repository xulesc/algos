package com.as.pms510;

import java.util.Comparator;

class DataKDNode extends KComparator<DataKDNode> {
	private float[] point_data;
	private float klass;

	DataKDNode(float[] p, float k) {
		point_data = p;
		klass = k;
	}
	
	public float[] getPointData() {
	    return point_data;
	}

	public float getKlass() {
		return klass;
	}

	@Override
	public String toString() {
		final StringBuffer sb = new StringBuffer();
		for (float f : point_data)
			sb.append(f + " ");
		return sb.toString();
	}

	@Override
	protected Comparator<DataKDNode> getComparator(final Integer axis) {
		return new Comparator<DataKDNode>() {
			@Override
			public int compare(DataKDNode arg0, DataKDNode arg1) {
				return Float.compare(arg0.point_data[axis],
						arg1.point_data[axis]);
			}
		};
	}

	@Override
	protected <T> Double sqDistance(T other) {
		final DataKDNode other_point = (DataKDNode) other;
		double sum = 0;
		for (int i = 0; i < point_data.length; i++)
			//if(point_data[i] != other_point.point_data[i])
			//	sum += 1;
			sum += Math.pow(point_data[i] - other_point.point_data[i], 2);
		return sum; //Math.pow(sum/point_data.length, 2);
	}

	@Override
	protected <T> Double axisSqDistance(T other, Integer axis) {
		final DataKDNode other_point = (DataKDNode) other;
		return Math.pow(point_data[axis] - other_point.point_data[axis], 2);
	}
}

