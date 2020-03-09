import dashing.Dashing;
import mash.Mash;

import java.io.IOException;
import java.util.Arrays;

/*
 *  Main.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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

public class Main {
    public static void main(String[] args) throws IOException {
        if (args.length==0){
            System.out.println("choose: mash, dashing, phylotest, or test (or both?)");
            System.exit(1);
        }
        switch (args[0]){
            case "mash":
                Mash.mash(Arrays.copyOfRange(args, 1, args.length));
                break;
            case "dashing":
                Dashing.dashing(Arrays.copyOfRange(args, 1, args.length));
                break;
            default:
                System.out.println("choose: mash, dashing, or test (or both?)");
                System.exit(1);
                break;
        }
    }
}
