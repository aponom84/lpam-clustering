package latna;

import Jama.Matrix;

import  smile.clustering.KMeans;
import java.util.HashMap;
import java.util.Vector;

public class KMeansWrapper implements EdgeClustering {
    @Override
    public HashMap<Integer, Vector<Integer>> clusterEdges(int k, Matrix distances, String suffix) {
        HashMap<Integer /*cluster*/, Vector<Integer> /*edges*/> clusters = new HashMap<>();

        KMeans km = new KMeans(distances.getArray(), k);
        final int[] results = km.getClusterLabel();

        for (int i = 0; i < distances.getColumnDimension(); ++i) {
            int cluster_index = results[i];
            clusters.computeIfAbsent(cluster_index, n -> new Vector<>());
            clusters.get(cluster_index).add(i + 1);
        }

        return clusters;
    }

    @Override
    public String getShortName() {
        return "kmn";
    }
}
