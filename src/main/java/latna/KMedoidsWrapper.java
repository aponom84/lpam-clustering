package latna;

import Jama.Matrix;

import smile.clustering.CLARANS;
import smile.math.distance.Distance;

import java.util.HashMap;
import java.util.Vector;

public class KMedoidsWrapper implements EdgeClustering{
    class MatrixDist implements Distance<Integer> {
        private Matrix dists;

        MatrixDist(Matrix dists)
        {
            this.dists = dists;
        }

        @Override
        public double d(Integer i, Integer j) {
            return this.dists.get(i, j);
        }
    }

    public HashMap<Integer /*cluster*/, Vector<Integer> /*edges*/> clusterEdges(int k, Matrix distances, String suffix) {
        HashMap<Integer /*cluster*/, Vector<Integer> /*edges*/> clusters = new HashMap<>();

        Vector<Integer> edges = new Vector<>();
        for (int i = 0; i < distances.getColumnDimension(); ++i) {
            edges.add(i);
        }

        CLARANS<Integer> km = new CLARANS<Integer>(edges.toArray(new Integer[0]), new MatrixDist(distances), k);
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
        return "KMD";
    }
}
