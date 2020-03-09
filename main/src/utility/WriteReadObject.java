package utility;

import dashing.DashingDistance;
import mash.MashDistance;

import java.io.*;
import java.util.ArrayList;

/*
 *  WriteReadObject.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 *
 * Benjamin Kaestle, 3.2020
 */

public class WriteReadObject {

    /**
     * writes a .phylip file for a distance matrix tableContent at a certain filepath
     * @param tableContent
     * @param filepath
     */
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


    /**
     * writes any seriazable object to a filepath
     * @param serObj
     * @param filepath
     */
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

    /**
     * reads any written object from filepath
     * @param filepath
     * @return
     */
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

    /**
     * writes a list of distances to a out filepath
     * @param distances
     * @param out
     */
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

    /**
     * writes a list of distances to a out filepath
     * @param distances
     * @param out
     */
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
