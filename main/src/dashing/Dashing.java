package dashing;

import utility.InputParameters;
import utility.WriteReadObject;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

public class Dashing {
    public static void main(String[] args){
        dashing(args);
    }

    private static void dashing(String[] args) {

        InputParameters parameters = new InputParameters();
        parameters.parseInputDashing(args);
        long startTime = System.currentTimeMillis();
        if (parameters.type.equals("sketch")){
            WriteReadObject.writeObjectToFile(dashingSketch(parameters),parameters.outputFile.concat(".dsk"));
        }
        else if (parameters.type.equals("dist")){
            DashingDistance[] dashingDistances = null;
            if (parameters.mashSketches !=null){
                dashingDistances = dashingDist(parameters.dashingSketches,parameters);
            }
            else{
                dashingDistances = dashingDist(dashingSketch(parameters),parameters);
            }
            if (parameters.tableOutput){
                printTable(tableOutput(dashingDistances,parameters.sequenceFiles));
            }else {
                DashingDistance[] distances1 = new DashingDistance[dashingDistances.length*2+parameters.sequenceFiles.size()];
                int i=0;
                for (i=0; i< dashingDistances.length; i++) {
                    distances1[2*i]= dashingDistances[i];
                    distances1[2*i+1]=new DashingDistance(dashingDistances[i]);
                }
                i*=2;
                for (String s : parameters.sequenceFiles){
                    distances1[i] = new DashingDistance(s,s,1,0,0,s,s,parameters.sketchSize);
                    i++;
                }
                Arrays.sort(distances1, Comparator.comparing(a -> a.getFilePath1()));
                Arrays.sort(distances1,Comparator.comparing(a -> a.getFilePath2()));
                System.out.println();
//                WriteReadObject.writeTxtFile(dashingDistances, parameters.outputFile);
                for (DashingDistance d : distances1) {

//                    System.out.print(d.getSameHashes()+", ");
                    d.printShort();
                }
                System.out.println();
            }
            if(parameters.outputPhylip){
                WriteReadObject.writePhylipFile(tableOutput(dashingDistances,parameters.sequenceFiles),parameters.outputFile);
            }
        }
        else if (parameters.type.equals("info")){
            for (DashingSketch dashingSketch : parameters.dashingSketches){
                System.out.println(dashingSketch.toString()+"\n");
            }
        }



//        mash_dist(mash_sketch(parameters),parameters)[0].print();
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        System.out.println(duration);

    }

    private static DashingDistance[] dashingDist(DashingSketch[] dashingSketches, InputParameters parameters) {
        return null;
    }

    private static DashingSketch[] dashingSketch(InputParameters parameters) {
        return null;
    }

    private static void printTable(String[][] tableOutput) {
    }

    private static String[][] tableOutput(DashingDistance[] dashingDistances, ArrayList<String> sequenceFiles) {
        return null;
    }
}
