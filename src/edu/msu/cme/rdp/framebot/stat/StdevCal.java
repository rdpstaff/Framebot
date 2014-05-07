package edu.msu.cme.rdp.framebot.stat;


import java.util.ArrayList;

/**
 *
 * @author wangqion
 */
public class StdevCal {
    public static class Std{
        private double totalCount;
        private double mean;
        private double stdev;
        Std(double t, double m, double s){
            totalCount = t;
            mean = m;
            stdev = s;
        }
        
        public double getTotalCount(){
            return totalCount;
        }
        public double getMean(){
            return mean;
        }
        public double getStdev(){
            return stdev;
        }

    }

    public static Std calStd(ArrayList<Double> valList){
        double totalCount = 0.0;
        double valSum = 0;
        for ( Double d: valList){
            if ( !d.isNaN()){
                valSum += d.doubleValue();
                totalCount++;
            }
        }
        double mean = valSum/totalCount;
        double sum = 0.0;
        for ( Double d: valList){
            if ( !d.isNaN()){
                sum += Math.pow(( d.doubleValue() - mean), 2);
            }
        }

        double stdev = Math.sqrt(sum/((double) totalCount -1));
        Std stdResult = new Std(totalCount, mean, stdev);
        return stdResult;
    }
    

    public static void main(String[] args){
        ArrayList<Double> valList = new ArrayList<Double>();
        valList.add(5.0);
        valList.add(12.0);
        valList.add(16.0);
        valList.add(21.0);
        valList.add(14.0);
        valList.add(15.0);
        Std result = StdevCal.calStd(valList);
        System.err.println("totalCount= " + result.totalCount + " mean= " + result.mean + " std=" + result.stdev);
    }

}
