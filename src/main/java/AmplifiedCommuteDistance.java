import Jama.Matrix;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.Graph;

import java.util.HashMap;

public class AmplifiedCommuteDistance implements DistanceMatrixCalculator {
    public Matrix calculateDistanceMatrix(Graph graph) {
        double a[][] = Utils.buildLineGraph(graph);

        int n = graph.getEdgeCount();

        double volG = 0;
        Matrix D = new Matrix(n, n);
        for (int i=0; i < n; i++) {
            double degree = 0;
            for (int j=0; j < n; j++  ) {
                degree = degree+a[i][j];
                volG = volG + a[i][j];
            }
            D.set(i, i, degree);
        }

        Matrix A = new Matrix(a);
        Matrix L = D.minus(A);

        L = Utils.pinv(L);

        Matrix d = new Matrix(n, n);

        for (int i = 0; i < n; i++) {
            d.set(i, i, 0);

            Matrix ei = new Matrix(1,n);
            ei.set(0, i, 1);
            for (int j = i+1; j < n; j++) {
                Matrix ej = new Matrix(1,n);
                ej.set(0, j, 1);

                Matrix left = ei.minus(ej);

                //left = left.times(volG);
                left =  left.times(L);

                Matrix right = ei.transpose().minus(ej.transpose());

                left = left.times(right);

                double R_ij = left.get(0, 0);

                double d_i = D.get(i, i);
                double d_j = D.get(j, j);
                double S_ij = R_ij - 1/d_i - 1/d_j;

                double u_ij = 2*a[i][j] / (d_i*d_j) - a[i][i] / (d_i*d_i) - a[j][j] / (d_j*d_j);
                double Camp_ij = Math.sqrt(S_ij + u_ij);

                d.set(i, j, Camp_ij);
                d.set(j, i, Camp_ij);
            }
        }
        return d;
    }

    @Override
    public String getShortName() {
        return "ACM";
    }
}
