package utility;

import java.io.*;

public class WriteReadObject {

    public static void writeObjectToFile(Object serObj, String filepath) {
        try {
            FileOutputStream fileOut = new FileOutputStream(filepath);
            ObjectOutputStream objectOut = new ObjectOutputStream(fileOut);
            objectOut.writeObject(serObj);
            objectOut.close();
            System.out.println("output file: "+filepath);

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
