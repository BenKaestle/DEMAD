import dashing.Dashing;
import mash.Mash;
import test.Testing;

import java.io.IOException;
import java.util.Arrays;

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
            case "phylotest":
                phylotree.Main.main(Arrays.copyOfRange(args, 1, args.length));
                break;
            case "test":
                Testing.testing(Arrays.copyOfRange(args, 1, args.length));
                break;
            default:
                System.out.println("choose: mash, dashing, or test (or both?)");
                System.exit(1);
                break;
        }
    }
}
