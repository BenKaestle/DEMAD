package utility;

import java.util.ArrayList;

/*
 *  SynchronizedSketch.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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

public class SynchronizedList<T> {
    ArrayList<T> list;

    public SynchronizedList(ArrayList<T> list) {
        this.list = list;
    }

    /**
     * get function for the list
     * @return
     */
    public synchronized T get(){
        if (list.isEmpty()) return null;
        return list.remove(0);
    }

    /**
     * returns weather the list is empty
     * @return
     */
    public synchronized boolean isEmpty(){
        return list.isEmpty();
    }

    /**
     * returns the current size of the list
     * @return
     */
    public synchronized int size(){
        return list.size();
    }
}
