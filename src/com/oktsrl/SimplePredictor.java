package com.oktsrl;

import java.util.HashMap;
import java.util.Map;

public class SimplePredictor {

	class AvgPredict{
		
		private double sum;
		int count;
		
		public AvgPredict() {
			sum=0d;
			count=0;
		}
		
		public void addContribution(double rank){
			sum+=rank;
			count++;
		}
		
		public double predict(){
			return sum / count;
		}
	}

	private Map<Object, AvgPredict> internalMap;
	
	public SimplePredictor() {
		this.internalMap=new HashMap<Object, AvgPredict>();
	}
	
	public void addContribution(Object id, double rank){
		if(internalMap.containsKey(id))
			internalMap.get(id).addContribution(rank);
		else{
			AvgPredict p = new AvgPredict();
			p.addContribution(rank);
			internalMap.put(id, p);
		}
	}
	
	public double predictRank(Object id){
		return this.internalMap.get(id).predict();
	}
	
}
