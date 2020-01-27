package utility;

import dashing.DashingDistance;
import mash.MashDistance;

import java.io.*;
import java.util.ArrayList;

public class WriteReadObject {

    public static int[] readJaccardIndex(String filepath){
        BufferedReader reader;
        int[] result = null;
        try {
            reader = new BufferedReader(new FileReader(filepath));
            String line = reader.readLine();
            int size = 0;
            while (line != null) {
                size++;
                line = reader.readLine();
            }
            result = new int[size];
            reader = new BufferedReader(new FileReader(filepath));
            line = reader.readLine();
            int i =0;
            while (line != null) {
                result[i] = Integer.parseInt(line.split("\\s+|\t")[line.split("\\s+|\t").length - 1].split("/")[0]);
                i++;
                line = reader.readLine();
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return result;
    }

    public static void writePhylipFile(String[][] tableContent, String filepath) {
        FileWriter fileWriter = null;
        String content = String.valueOf(tableContent.length - 1);
        int rowNumber = 0;
        boolean firstInRow;
        boolean firstRow = true;
        String blanks;
        String title;
        String titleList = "";
        for (String[] row : tableContent) {
            firstInRow = true;
            for (String cell : row) {
                if (firstRow) {

                } else if (firstInRow) {
                    firstInRow = false;
                    title = "G" + rowNumber;
                    titleList += title + "\t" + cell + "\n";
                    blanks = new String(new char[10 - title.length()]).replace("\0", " ");
                    content += "\n" + title + blanks;
                } else content += "  " + cell;
            }
            rowNumber++;
            firstRow = false;
        }
        try {
            fileWriter = new FileWriter(filepath.concat(".phylip"));
            fileWriter.write(content);
            fileWriter.close();
            fileWriter = new FileWriter(filepath.concat("_phylip_help.txt"));
            fileWriter.write(titleList);
            fileWriter.close();
            System.out.println("output files: " + filepath.concat(".phylip") + ", " + filepath.concat("_phylip_help.txt"));
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static String[][] readTSVtable(String table) {
        int columns = 0;
        int rows = 0;
        for (String row : table.split("\n")) {
            if (rows == 0) {
                for (String cell : row.split("\\s+|\t")) {
                    columns++;
                }
            }
            rows++;
        }
        String[][] result = new String[rows][columns];
        columns = 0;
        rows = 0;
        for (String row : table.split("\n")) {
            columns=0;
            for (String cell : row.split("\\s+|\t")) {
                result[rows][columns] = cell;
                columns++;
            }
            rows++;
        }
        return result;
    }

    public static void writeObjectToFile(Object serObj, String filepath) {
        try {
            FileOutputStream fileOut = new FileOutputStream(filepath);
            ObjectOutputStream objectOut = new ObjectOutputStream(fileOut);
            objectOut.writeObject(serObj);
            objectOut.close();
            System.out.println("output file: " + filepath);

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static Object readObjectFromFile(String filepath) {
        try {
            FileInputStream fin = new FileInputStream(filepath);
            ObjectInputStream ois = new ObjectInputStream(fin);
            return ois.readObject();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static ArrayList<String> readTxt(String filepath){
        BufferedReader reader;
        ArrayList<String> result = new ArrayList<>();
        try {
            reader = new BufferedReader(new FileReader(filepath));
            String line = reader.readLine();
            while (line != null) {
                result.add(line);
                line = reader.readLine();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return result;
    }

    public static void writeTxt (String filepath, ArrayList<String> data){
        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(filepath+".txt");
            for (String d : data){
                fileWriter.write(d+"\n");
            }
            fileWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void writeTxtFile(MashDistance[] distances, String out) {
        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(out+".txt");
            for (MashDistance d : distances){
                fileWriter.write(d.toString()+"\n");
            }
            fileWriter.close();
            System.out.println("output file: "+out+".txt");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void writeTxtFile(DashingDistance[] distances, String out) {
        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(out+".txt");
            for (DashingDistance d : distances){
                fileWriter.write(d.toString()+"\n");
            }
            fileWriter.close();
            System.out.println("output file: "+out+".txt");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
