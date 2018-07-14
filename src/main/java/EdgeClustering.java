import Jama.Matrix;

import java.util.HashMap;
import java.util.Vector;

public interface EdgeClustering {
    HashMap<Integer /*cluster*/, Vector<Integer> /*edges*/> clusterEdges(int k, Matrix distances, String suffix);
    String getShortName();
}
