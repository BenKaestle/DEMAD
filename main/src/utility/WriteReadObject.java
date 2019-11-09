package utility;

import java.io.*;

public class WriteReadObject {

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
                for (String cell : row.split("\t")) {
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
            for (String cell : row.split("\t")) {
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
}
